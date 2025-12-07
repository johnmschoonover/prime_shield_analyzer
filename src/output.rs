#![allow(clippy::manual_is_multiple_of)]
use crate::config::Config;
use crate::stats::Statistics;
use csv::Writer;
use serde::Serialize;
use std::collections::{BTreeMap, BTreeSet};
use std::error::Error;
use std::fs;
use std::path::Path;

pub fn write_results(
    stats: &Statistics,
    config: &Config,
    max_n: u64,
) -> Result<(), Box<dyn Error>> {
    fs::create_dir_all(&config.output_dir)?;

    write_global_stats(stats, config)?;
    write_gap_spectrum(stats, config, max_n)?;
    write_oscillation_series(stats, config)?;

    Ok(())
}

#[derive(Serialize)]
struct GlobalStatsRecord {
    total_primes_p: u64,
    total_primes_s: u64,
    global_ratio_s_p: f64,
}

fn write_global_stats(stats: &Statistics, config: &Config) -> Result<(), Box<dyn Error>> {
    let path = Path::new(&config.output_dir).join("global_stats.csv");
    let mut wtr = Writer::from_path(path)?;

    let ratio = if stats.total_primes > 0 {
        stats.total_s_primes as f64 / stats.total_primes as f64
    } else {
        0.0
    };

    let record = GlobalStatsRecord {
        total_primes_p: stats.total_primes,
        total_primes_s: stats.total_s_primes,
        global_ratio_s_p: ratio,
    };

    wtr.serialize(record)?;
    wtr.flush()?;
    Ok(())
}

#[derive(Debug)]
struct ShieldingInfo {
    shield_score: u32,
    shield_primes: String,
    theoretical_boost: f64,
}

// Helper to get unique prime factors
fn get_prime_factors(mut n: u64) -> Vec<u64> {
    let mut factors = Vec::new();
    if n < 2 {
        return factors;
    }

    // Handle 2 separately
    if n % 2 == 0 {
        factors.push(2);
        while n % 2 == 0 {
            n /= 2;
        }
    }

    // Handle odd factors
    let mut i = 3;
    while i * i <= n {
        if n % i == 0 {
            factors.push(i);
            while n % i == 0 {
                n /= i;
            }
        }
        i += 2;
    }
    if n > 1 {
        factors.push(n);
    }
    factors
}

fn calculate_shielding_info(g: u64) -> ShieldingInfo {
    let mut unique_shields = BTreeSet::new();

    // 1. Neighbor Hazards (g - 1)
    // Corresponds to g = 1 mod q (Natural Shield)
    for p in get_prime_factors(g - 1) {
        unique_shields.insert(p);
    }

    // 2. Neighbor Hazards (g + 1)
    // Corresponds to g = -1 mod q (Selection Shield)
    for p in get_prime_factors(g + 1) {
        unique_shields.insert(p);
    }

    // Filter out 2 (handled by parity)
    unique_shields.remove(&2);

    let mut theoretical_boost = 1.0;
    let mut shield_primes_vec: Vec<u64> = Vec::new();

    for &q in &unique_shields {
        shield_primes_vec.push(q);
        theoretical_boost *= q as f64 / (q as f64 - 1.0);
    }

    // Mod 3 Trap: If g % 3 == 0, S = 2p + g - 1 fails whenever p = 2 mod 3 (50% of odd primes).
    // Baseline probability of a random odd number being coprime to 3 is 2/3 (66%).
    // The ratio of Trap Success (1/2) to Baseline (2/3) is (1/2) / (2/3) = 3/4 = 0.75.
    if g % 3 == 0 {
        theoretical_boost *= 0.75;
    }

    let shield_primes = shield_primes_vec
        .iter()
        .map(|p| p.to_string())
        .collect::<Vec<String>>()
        .join(",");

    ShieldingInfo {
        shield_score: unique_shields.len() as u32,
        shield_primes,
        theoretical_boost,
    }
}

#[derive(Serialize)]
struct GapSpectrumRecord {
    gap_size: u64,
    count: u64,
    successes: u64,
    success_rate: f64,
    expected_rate_heuristic: f64,
    shield_score: u32,
    shield_primes: String,
    theoretical_boost: f64,
}

fn write_gap_spectrum(
    stats: &Statistics,
    config: &Config,
    max_n: u64,
) -> Result<(), Box<dyn Error>> {
    let path = Path::new(&config.output_dir).join("gap_spectrum.csv");
    let mut wtr = Writer::from_path(path)?;

    // Updated Heuristic: 2.0 / ln(N) because we scan only odd numbers
    // This provides a more accurate baseline for prime density in this context.
    let expected_rate = 2.0 / (max_n as f64).ln();

    let sorted_gaps: BTreeMap<_, _> = stats.gap_spectrum.iter().collect();

    for (&gap_size, &(count, successes)) in sorted_gaps {
        let success_rate = if count > 0 {
            successes as f64 / count as f64
        } else {
            0.0
        };
        let shielding_info = calculate_shielding_info(gap_size);

        let record = GapSpectrumRecord {
            gap_size,
            count,
            successes,
            success_rate,
            expected_rate_heuristic: expected_rate,
            shield_score: shielding_info.shield_score,
            shield_primes: shielding_info.shield_primes,
            theoretical_boost: shielding_info.theoretical_boost,
        };
        wtr.serialize(record)?;
    }

    wtr.flush()?;
    Ok(())
}

fn write_oscillation_series(stats: &Statistics, config: &Config) -> Result<(), Box<dyn Error>> {
    let path = Path::new(&config.output_dir).join("oscillation_series.csv");
    let mut wtr = Writer::from_path(path)?;

    // Dynamically build headers
    let mut headers: Vec<String> = vec![
        "bin_start".to_string(),
        "bin_end".to_string(),
        "prime_count_p".to_string(),
        "prime_count_s".to_string(),
        "ratio_s_p".to_string(),
    ];
    for &g in &stats.target_gaps {
        headers.push(format!("gap_{}_rate", g));
    }
    wtr.write_record(headers.iter())?;

    // Dynamically write data rows
    for bin in &stats.bins {
        // Only write bins that contain actual prime data (prime_count_p > 0)
        if bin.prime_count_p == 0 {
            continue;
        }

        let ratio_s_p = if bin.prime_count_p > 0 {
            bin.prime_count_s as f64 / bin.prime_count_p as f64
        } else {
            0.0
        };

        let mut record: Vec<String> = vec![
            bin.bin_start.to_string(),
            bin.bin_end.to_string(),
            bin.prime_count_p.to_string(),
            bin.prime_count_s.to_string(),
            ratio_s_p.to_string(),
        ];

        for &g in &stats.target_gaps {
            let occurrences = bin.gap_occurrences.get(&g).cloned().unwrap_or(0);
            let successes = bin.gap_successes.get(&g).cloned().unwrap_or(0);
            let rate = if occurrences > 0 {
                successes as f64 / occurrences as f64
            } else {
                0.0
            };
            record.push(rate.to_string());
        }
        wtr.write_record(record.iter())?;
    }
    wtr.flush()?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_shielding_logic() {
        // Gap 2
        // Factors(2) -> 2 (Filtered)
        // Factors(1) -> []
        // Factors(3) -> 3
        // Result: 3. Boost 1.5.
        let info_2 = calculate_shielding_info(2);
        assert_eq!(info_2.shield_score, 1);
        assert_eq!(info_2.shield_primes, "3");
        assert_eq!(info_2.theoretical_boost, 1.5);

        // Gap 4
        // Factors(4) -> 2 (Filtered)
        // Factors(3) -> 3
        // Factors(5) -> 5
        // Result: 3, 5. Boost 1.5 * 1.25
        let info_4 = calculate_shielding_info(4);
        assert_eq!(info_4.shield_score, 2);
        assert_eq!(info_4.shield_primes, "3,5");
        assert_eq!(info_4.theoretical_boost, 1.5 * 1.25);

        // Gap 6
        // Factors(6) -> 2, 3 (Wheel - IGNORED)
        // Factors(5) -> 5
        // Factors(7) -> 7
        // Result: 5, 7. Boost 1.25 * 1.166 * 0.75 (Mod 3 Penalty)
        let info_6 = calculate_shielding_info(6);
        assert_eq!(info_6.shield_score, 2);
        assert_eq!(info_6.shield_primes, "5,7");
        assert_eq!(info_6.theoretical_boost, 1.25 * (7.0 / 6.0) * 0.75);

        // Gap 30
        // Factors(30) -> 2, 3, 5 (Wheel - IGNORED)
        // Factors(29) -> 29
        // Factors(31) -> 31
        // Result: 29, 31. Boost * 0.75 (Mod 3 Penalty)
        let info_30 = calculate_shielding_info(30);
        assert_eq!(info_30.shield_primes, "29,31");
        assert_eq!(info_30.shield_score, 2);
        assert_eq!(
            info_30.theoretical_boost,
            (29.0 / 28.0) * (31.0 / 30.0) * 0.75
        );
    }
}
