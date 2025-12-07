#![allow(clippy::manual_is_multiple_of)]
use crate::config::Config;
use crate::stats::Statistics;
use csv::Writer;
use serde::Serialize;
use std::collections::BTreeMap;
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

// Pre-compute primes up to 100 for the shielding calculation.
const SMALL_PRIMES: &[u32] = &[
    3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97,
];

fn calculate_shielding_info(g: u64) -> ShieldingInfo {
    // --- DEBUG TRAP START ---
    if g == 56 {
        println!("!!! DEBUG TRAP FOR GAP 56 !!!");
        println!("Mod 3 Residue: {} (Should be 2)", g % 3);
        println!("Mod 19 Residue: {} (Should be 18)", g % 19);
    }
    // --- DEBUG TRAP END ---

    let mut shield_score = 0;
    let mut shield_primes_vec = Vec::new();
    let mut theoretical_boost = 1.0;

    // The Mod 3 Rule
    // If g % 3 == 1, S = 2p is never 0 mod 3.
    // If g % 3 == 2, valid prime gaps only exist for p = 2 mod 3.
    // (p = 1 mod 3 implies p+g is divisible by 3, so neighbor is composite).
    // In this forced case, S = 2(2) + 2 - 1 = 5 = 2 mod 3.
    // Thus, ANY gap not divisible by 3 is shielded from 3.
    if g % 3 != 0 {
        shield_score += 1;
        shield_primes_vec.push(3);
        theoretical_boost *= 3.0 / 2.0;
    }

    // The General Rule (q >= 5)
    for &q in SMALL_PRIMES.iter().skip(1) {
        // Skip 3 as it's already handled
        let q_u64 = q as u64;
        // CORRECTION: A gap is shielded if g = 1 mod q.
        // Proof: S = 2p + g - 1. If g = 1 mod q, then g - 1 = k*q.
        // S = 2p + k*q = 2p (mod q). Since p is prime > q, 2p != 0 (mod q).
        // Thus S is never divisible by q.
        if g % q_u64 == 1 {
            shield_score += 1;
            shield_primes_vec.push(q);
            theoretical_boost *= q_u64 as f64 / (q_u64 - 1) as f64;
        }
    }

    let shield_primes = shield_primes_vec
        .iter()
        .map(|p| p.to_string())
        .collect::<Vec<String>>()
        .join(",");

    ShieldingInfo {
        shield_score,
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

    let expected_rate = 1.0 / (max_n as f64).ln();

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
        // This prevents long stretches of zero data at the end of the graph
        // for smaller max_exponent values, improving graph readability.
        if bin.prime_count_p == 0 {
            continue;
        }

        let ratio_s_p = bin.prime_count_s as f64 / bin.prime_count_p as f64;

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
        // Test Gap 2: Shielded by 3 (Selection Bias)
        let info_2 = calculate_shielding_info(2);
        assert_eq!(info_2.shield_score, 1);
        assert_eq!(info_2.shield_primes, "3");
        assert_eq!(info_2.theoretical_boost, 1.5);

        // Test Gap 4: Shielded by 3 only (4 % 5 = 4 != 1)
        let info_4 = calculate_shielding_info(4);
        assert_eq!(info_4.shield_score, 1);
        assert_eq!(info_4.shield_primes, "3");
        assert_eq!(info_4.theoretical_boost, 1.5);

        // Test Gap 6: Shielded by 5 (6 % 5 = 1)
        let info_6 = calculate_shielding_info(6);
        assert_eq!(info_6.shield_score, 1);
        assert_eq!(info_6.shield_primes, "5");
        assert_eq!(info_6.theoretical_boost, 1.25);

        // Test Gap 34: Shielded by 3 and 11
        let info_34 = calculate_shielding_info(34);
        assert_eq!(info_34.shield_score, 2);
        assert_eq!(info_34.shield_primes, "3,11");
        assert_eq!(info_34.theoretical_boost, (3.0 / 2.0) * (11.0 / 10.0));

        // Test Gap 56: Shielded by 3, 5, 11
        let info_56 = calculate_shielding_info(56);
        assert_eq!(info_56.shield_score, 3);
        assert_eq!(info_56.shield_primes, "3,5,11");
        assert_eq!(
            info_56.theoretical_boost,
            (3.0 / 2.0) * (5.0 / 4.0) * (11.0 / 10.0)
        );
    }
}
