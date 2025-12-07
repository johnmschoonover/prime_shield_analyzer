mod config;
mod output;
mod report;
mod sieve;
mod stats;

use clap::Parser;
use indicatif::ProgressBar;
use sieve::{PrimalityChecker, PrimeIterator};
use stats::Statistics;

use crate::config::Config;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let config = Config::parse();

    // --- Config Validation ---
    if config.gaps.is_empty() {
        eprintln!("Error: No gap sizes provided.");
        std::process::exit(1);
    }
    for &gap in &config.gaps {
        if gap == 0 {
            eprintln!("Error: Gap size cannot be 0.");
            std::process::exit(1);
        }
        if gap % 2 != 0 && gap != 1 {
            eprintln!("Error: Gap size {} is odd.", gap);
            std::process::exit(1);
        }
    }

    // --- Optimization 1: Lookup Table for Target Gaps ---
    // Instead of a HashSet (Hashing overhead), use a boolean array for small gaps.
    // Most target gaps are < 256.
    let max_target = *config.gaps.iter().max().unwrap_or(&0) as usize;
    let use_lookup = max_target < 1024; // Arbitrary limit for array size
    let mut target_lookup = if use_lookup {
        vec![false; max_target + 1]
    } else {
        Vec::new()
    };
    let target_set: std::collections::HashSet<u64> = if !use_lookup {
        config.gaps.iter().cloned().collect()
    } else {
        for &g in &config.gaps {
            target_lookup[g as usize] = true;
        }
        std::collections::HashSet::new()
    };

    // Quick closure to check if a gap is "interesting" without hashing
    let is_target_gap = |g: u64| -> bool {
        if use_lookup {
            if g as usize <= max_target {
                target_lookup[g as usize]
            } else {
                false
            }
        } else {
            target_set.contains(&g)
        }
    };

    let mut sorted_target_gaps = config.gaps.clone();
    sorted_target_gaps.sort_unstable();

    let max_n = 10u64.pow(config.max_exponent);
    let segment_size_bytes = config.segment_size_kb * 1024;

    println!("Max N (10^{}): {}", config.max_exponent, max_n);
    println!("Bins: {}", config.bins);
    println!("Output Dir: {}", config.output_dir);
    println!("Tracking Gaps: {:?}", sorted_target_gaps);

    let mut prime_iterator = PrimeIterator::new(max_n, segment_size_bytes);
    let analysis_limit = max_n * 2;
    let mut primality_checker = PrimalityChecker::new(analysis_limit, segment_size_bytes);

    let mut stats = Statistics::new(max_n, config.bins, &sorted_target_gaps);

    let bar = ProgressBar::new(max_n);
    bar.set_style(indicatif::ProgressStyle::default_bar()
        .template("{spinner:.green} [{elapsed_precise}] [{bar:40.cyan/blue}] {pos:>7}/{len:7} ({eta}) {msg}")?
        .progress_chars("#>-"));

    let mut p_prev = 2;

    // Handle first prime (2) manually
    if max_n >= 2 {
        stats.total_primes += 1;
        if let Some(bin_index) = stats.get_bin_index(2) {
            stats.bins[bin_index].prime_count_p += 1;
        }
        prime_iterator.next();
    }

    // --- Optimization 2: Fast-Path Gap Cache ---
    // Avoid HashMap insertions for every single gap. Buffer small gaps in an array.
    const MAX_FAST_GAP: usize = 320;
    let mut gap_counts = vec![0u64; MAX_FAST_GAP];
    let mut gap_successes = vec![0u64; MAX_FAST_GAP];

    // Main Loop
    for p_current in prime_iterator {
        stats.total_primes += 1;

        // Hot Path: Bin Index Calculation
        // Inlining this logic manually or trusting the compiler.
        // For simple bins, (n / size) is fast.
        if let Some(bin_index) = stats.get_bin_index(p_current) {
            stats.bins[bin_index].prime_count_p += 1;

            let gap = p_current - p_prev;
            let s = p_current + p_prev - 1;

            // 1. Record Gap Occurrence
            if (gap as usize) < MAX_FAST_GAP {
                gap_counts[gap as usize] += 1;
            } else {
                stats.gap_spectrum.entry(gap).or_insert((0, 0)).0 += 1;
            }

            // High-interest gap tracking (Bin specific)
            if is_target_gap(gap) {
                *stats.bins[bin_index]
                    .gap_occurrences
                    .entry(gap)
                    .or_insert(0) += 1;
            }

            // 2. Check S
            if primality_checker.is_prime(s) {
                stats.total_s_primes += 1;
                stats.bins[bin_index].prime_count_s += 1;

                if (gap as usize) < MAX_FAST_GAP {
                    gap_successes[gap as usize] += 1;
                } else {
                    stats.gap_spectrum.entry(gap).or_insert((0, 0)).1 += 1;
                }

                if is_target_gap(gap) {
                    *stats.bins[bin_index].gap_successes.entry(gap).or_insert(0) += 1;
                }
            }
        }

        p_prev = p_current;

        // --- Optimization 3: Throttled UI Updates ---
        // Only update progress bar every ~250k primes.
        // 0x3FFFF is a bitmask for 262,143.
        if (stats.total_primes & 0x3FFFF) == 0 {
            bar.set_position(p_current);
        }
    }

    bar.finish_with_message("Analysis complete.");

    // Flush Fast-Path Cache to HashMap
    for gap in 0..MAX_FAST_GAP {
        if gap_counts[gap] > 0 {
            let entry = stats.gap_spectrum.entry(gap as u64).or_insert((0, 0));
            entry.0 += gap_counts[gap];
            entry.1 += gap_successes[gap];
        }
    }

    println!("Writing results...");
    output::write_results(&stats, &config, max_n)?;

    if config.web_report {
        println!("Generating HTML report...");
        report::generate_report(&config, max_n)?;
        println!("Report generated at {}/index.html", config.output_dir);
    }

    Ok(())
}
