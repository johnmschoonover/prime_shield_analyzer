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

    // Validate gaps
    if config.gaps.is_empty() {
        eprintln!("Error: No gap sizes provided. Please provide at least one gap size.");
        std::process::exit(1);
    }
    for &gap in &config.gaps {
        if gap == 0 {
            eprintln!("Error: Gap size cannot be 0.");
            std::process::exit(1);
        }
        if gap % 2 != 0 && gap != 1 {
            // Allow gap of 1 only for special cases, but statistical analysis expects even.
            eprintln!(
                "Error: Gap size {} is odd. All relevant prime gaps (except for the first, between 2 and 3) are even. Please provide even gap sizes.",
                gap
            );
            std::process::exit(1);
        }
    }

    let target_gaps_set: std::collections::HashSet<u64> = config.gaps.iter().cloned().collect();
    let mut sorted_target_gaps = config.gaps.clone();
    sorted_target_gaps.sort_unstable(); // For consistent CSV/HTML output order

    let max_n = 10u64.pow(config.max_exponent);

    // Use the user-defined segment size, converting from KB to Bytes.
    let segment_size_bytes = config.segment_size_kb * 1024;

    println!("Max N (10^{}): {}", config.max_exponent, max_n);
    println!("Bins: {}", config.bins);
    println!("Output Dir: {}", config.output_dir);
    println!("Using Segment Size: {} KB", config.segment_size_kb);
    println!("Tracking Gaps: {:?}", sorted_target_gaps);

    // The sieve for generating p_n only needs to go up to max_n.
    let mut prime_iterator = PrimeIterator::new(max_n, segment_size_bytes);

    // The checker needs to handle sums S = p_n + p_{n+1} - 1.
    // So S can be close to 2 * max_n.
    let analysis_limit = max_n * 2;
    let mut primality_checker = PrimalityChecker::new(analysis_limit, segment_size_bytes);

    let mut stats = Statistics::new(max_n, config.bins, &sorted_target_gaps);

    let bar = ProgressBar::new(max_n);
    bar.set_style(indicatif::ProgressStyle::default_bar()
        .template("{spinner:.green} [{elapsed_precise}] [{bar:40.cyan/blue}] {pos:>7}/{len:7} ({eta})")?
        .progress_chars("#>-"));

    let mut p_prev = 2; // The first prime

    // Manually handle the first prime (2) since our loop starts with the second one
    if max_n >= 2 {
        stats.total_primes += 1;
        if let Some(bin_index) = stats.get_bin_index(2) {
            stats.bins[bin_index].prime_count_p += 1;
        }
        prime_iterator.next(); // Consume '2' from iterator
    }

    for p_current in prime_iterator {
        stats.total_primes += 1;

        // Update stats for p_current
        if let Some(bin_index) = stats.get_bin_index(p_current) {
            stats.bins[bin_index].prime_count_p += 1;
        }

        // Calculate gap and S
        let gap = p_current - p_prev;
        let s = p_current + p_prev - 1;

        // Update gap spectrum (occurrences)
        stats.gap_spectrum.entry(gap).or_insert((0, 0)).0 += 1;

        // Update high-interest gap occurrences in the correct bin
        // The occurrence is tied to the location of p_current
        if target_gaps_set.contains(&gap) {
            if let Some(bin_index) = stats.get_bin_index(p_current) {
                *stats.bins[bin_index]
                    .gap_occurrences
                    .entry(gap)
                    .or_insert(0) += 1;
            }
        }

        // Check if S is prime
        if primality_checker.is_prime(s) {
            stats.total_s_primes += 1;

            // Update gap spectrum (successes)
            stats.gap_spectrum.entry(gap).or_insert((0, 0)).1 += 1;

            // Update bin stats for S
            if let Some(bin_index) = stats.get_bin_index(s) {
                stats.bins[bin_index].prime_count_s += 1;
            }

            // Update high-interest gap successes in the correct bin
            // The success is also tied to the location of p_current
            if target_gaps_set.contains(&gap) {
                if let Some(bin_index) = stats.get_bin_index(p_current) {
                    *stats.bins[bin_index].gap_successes.entry(gap).or_insert(0) += 1;
                }
            }
        }

        p_prev = p_current;
        bar.set_position(p_current);
    }
    bar.finish_with_message("Sieving and analysis complete.");

    println!("Writing results to disk...");
    output::write_results(&stats, &config, max_n)?;
    println!("Done.");

    if config.web_report {
        println!("Generating HTML report...");
        report::generate_report(&config, max_n)?;
        println!("Report generated at {}/index.html", config.output_dir);
    }

    Ok(())
}
