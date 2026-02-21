mod config;
mod output;
mod report;
mod sieve;
mod stats;

use clap::Parser;
use indicatif::ProgressBar;
use rayon::prelude::*;
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
    let primality_checker = PrimalityChecker::new(analysis_limit, segment_size_bytes);

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

    // --- Parallel Batch Processing ---
    const BATCH_SIZE: usize = 262_144; // 2^18
    let mut batch = Vec::with_capacity(BATCH_SIZE);

    for p_current in prime_iterator {
        batch.push(p_current);

        if batch.len() >= BATCH_SIZE {
            process_batch(
                &batch,
                p_prev,
                &primality_checker,
                &mut stats,
                max_n,
                config.bins,
                &sorted_target_gaps,
                &is_target_gap,
            );

            p_prev = *batch.last().unwrap();
            bar.inc(batch.len() as u64);
            batch.clear();
        }
    }

    // Process remaining primes
    if !batch.is_empty() {
        process_batch(
            &batch,
            p_prev,
            &primality_checker,
            &mut stats,
            max_n,
            config.bins,
            &sorted_target_gaps,
            &is_target_gap,
        );
        bar.inc(batch.len() as u64);
    }

    bar.finish_with_message("Analysis complete.");

    println!("Writing results...");
    output::write_results(&stats, &config, max_n)?;

    if config.web_report {
        println!("Generating HTML report...");
        report::generate_report(&config, max_n)?;
        println!("Report generated at {}/index.html", config.output_dir);
    }

    Ok(())
}

#[allow(clippy::too_many_arguments)]
fn process_batch<F>(
    batch: &[u64],
    p_prev_start: u64,
    primality_checker: &PrimalityChecker,
    global_stats: &mut Statistics,
    max_n: u64,
    num_bins: usize,
    target_gaps: &[u64],
    is_target_gap: &F,
) where
    F: Fn(u64) -> bool + Sync + Send,
{
    // Parallel Map-Reduce
    let batch_stats = batch
        .par_iter()
        .zip(rayon::iter::once(&p_prev_start).chain(batch.par_iter()))
        .fold(
            || Statistics::new(max_n, num_bins, target_gaps),
            |mut local_stats, (&p_current, &p_prev)| {
                local_stats.total_primes += 1;

                if let Some(bin_index) = local_stats.get_bin_index(p_current) {
                    local_stats.bins[bin_index].prime_count_p += 1;

                    let gap = p_current - p_prev;
                    let s = p_current + p_prev - 1;

                    // 1. Record Gap Occurrence
                    local_stats.gap_spectrum.entry(gap).or_insert((0, 0)).0 += 1;

                    if is_target_gap(gap) {
                        *local_stats.bins[bin_index]
                            .gap_occurrences
                            .entry(gap)
                            .or_insert(0) += 1;
                    }

                    // 2. Check S
                    if primality_checker.is_prime(s) {
                        local_stats.total_s_primes += 1;
                        local_stats.bins[bin_index].prime_count_s += 1;

                        local_stats.gap_spectrum.entry(gap).or_insert((0, 0)).1 += 1;

                        if is_target_gap(gap) {
                            *local_stats.bins[bin_index]
                                .gap_successes
                                .entry(gap)
                                .or_insert(0) += 1;
                        }
                    }
                }
                local_stats
            },
        )
        .reduce(
            || Statistics::new(max_n, num_bins, target_gaps),
            |mut a, b| {
                a.merge(&b);
                a
            },
        );

    global_stats.merge(&batch_stats);
}
