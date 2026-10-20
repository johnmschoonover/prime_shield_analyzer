use crate::config::Config;
use crate::stats::Statistics;
use csv::Writer;
use serde::Serialize;
use std::collections::BTreeMap;
use std::error::Error;
use std::fs;
use std::path::Path;
use std::sync::OnceLock;

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

/// Largest prime included in the local-factor product. The neutral factor is
/// 1 - 2/q^2 + O(1/q^3), so the tail beyond this bound changes the boost by
/// well under 0.1%.
const LOCAL_FACTOR_PRIME_BOUND: u64 = 1000;

fn odd_primes() -> &'static [u64] {
    static PRIMES: OnceLock<Vec<u64>> = OnceLock::new();
    PRIMES.get_or_init(|| {
        let limit = LOCAL_FACTOR_PRIME_BOUND as usize;
        let mut sieve = vec![true; limit + 1];
        sieve[0] = false;
        sieve[1] = false;
        let mut i = 2;
        while i * i <= limit {
            if sieve[i] {
                let mut j = i * i;
                while j <= limit {
                    sieve[j] = false;
                    j += i;
                }
            }
            i += 1;
        }
        (3..=limit)
            .filter(|&n| sieve[n])
            .map(|n| n as u64)
            .collect()
    })
}

/// How a single odd prime `q` interacts with the sum S = 2p + g - 1 when both
/// `p` and `p + g` are prime.
///
/// Both parents are prime, so `p` avoids the residues 0 and -g (mod q). The
/// sum is divisible by q exactly when p = (1 - g)/2 (mod q). Three cases:
///
/// * **Shield** (g = +1 or -1 mod q): the bad residue coincides with one the
///   parents already exclude, so S is never divisible by q. Relative to a
///   random odd integer, which survives q with probability (q-1)/q, the sum
///   survives with probability 1, a boost of q/(q-1).
/// * **Trap** (g = 0 mod q): the parents exclude only one residue instead of
///   two, so S is divisible by q for 1 of the q-1 remaining residues. The
///   survival ratio is ((q-2)/(q-1)) / ((q-1)/q) = q(q-2)/(q-1)^2, a penalty.
/// * **Neutral** (anything else): S is divisible by q for 1 of q-2 residues,
///   giving ((q-3)/(q-2)) / ((q-1)/q) = q(q-3)/((q-1)(q-2)), a mild penalty.
///
/// The product of these factors over all odd primes is the Bateman-Horn
/// local correction for the linear system (n, n+g, 2n+g-1), conditioned on
/// the first two forms being prime.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum LocalEffect {
    Shield,
    Trap,
    Neutral,
}

pub fn local_effect(g: u64, q: u64) -> LocalEffect {
    match g % q {
        0 => LocalEffect::Trap,
        r if r == 1 || r == q - 1 => LocalEffect::Shield,
        _ => LocalEffect::Neutral,
    }
}

pub fn local_factor(g: u64, q: u64) -> f64 {
    let qf = q as f64;
    match local_effect(g, q) {
        LocalEffect::Shield => qf / (qf - 1.0),
        LocalEffect::Trap => qf * (qf - 2.0) / ((qf - 1.0) * (qf - 1.0)),
        LocalEffect::Neutral => qf * (qf - 3.0) / ((qf - 1.0) * (qf - 2.0)),
    }
}

#[derive(Debug)]
pub struct ShieldingInfo {
    pub shield_primes: Vec<u64>,
    pub trap_primes: Vec<u64>,
    /// Predicted success rate divided by the size-matched baseline 2/ln(S).
    pub theoretical_boost: f64,
}

pub fn calculate_shielding_info(g: u64) -> ShieldingInfo {
    let mut shield_primes = Vec::new();
    let mut trap_primes = Vec::new();
    let mut theoretical_boost = 1.0;

    // The only odd gap is the (2, 3) pair. The residue argument assumes an
    // odd prime p, so report that single pair as unmodelled.
    if g % 2 == 1 {
        return ShieldingInfo {
            shield_primes,
            trap_primes,
            theoretical_boost,
        };
    }

    for &q in odd_primes() {
        match local_effect(g, q) {
            LocalEffect::Shield => shield_primes.push(q),
            LocalEffect::Trap => trap_primes.push(q),
            LocalEffect::Neutral => {}
        }
        theoretical_boost *= local_factor(g, q);
    }

    ShieldingInfo {
        shield_primes,
        trap_primes,
        theoretical_boost,
    }
}

fn join_primes(primes: &[u64]) -> String {
    primes
        .iter()
        .map(|p| p.to_string())
        .collect::<Vec<String>>()
        .join(",")
}

#[derive(Serialize)]
struct GapSpectrumRecord {
    gap_size: u64,
    count: u64,
    successes: u64,
    success_rate: f64,
    /// Mean of 2/ln(S) over the pairs with this gap: what a random odd
    /// integer of the same size would score.
    baseline_rate: f64,
    /// baseline_rate * theoretical_boost: the model's prediction.
    predicted_rate: f64,
    theoretical_boost: f64,
    shield_score: usize,
    shield_primes: String,
    trap_primes: String,
}

fn write_gap_spectrum(
    stats: &Statistics,
    config: &Config,
    _max_n: u64,
) -> Result<(), Box<dyn Error>> {
    let path = Path::new(&config.output_dir).join("gap_spectrum.csv");
    let mut wtr = Writer::from_path(path)?;

    let sorted_gaps: BTreeMap<_, _> = stats.gap_spectrum.iter().collect();

    for (&gap_size, tally) in sorted_gaps {
        let (success_rate, baseline_rate) = if tally.count > 0 {
            (
                tally.successes as f64 / tally.count as f64,
                tally.baseline_sum / tally.count as f64,
            )
        } else {
            (0.0, 0.0)
        };
        let info = calculate_shielding_info(gap_size);

        let record = GapSpectrumRecord {
            gap_size,
            count: tally.count,
            successes: tally.successes,
            success_rate,
            baseline_rate,
            predicted_rate: baseline_rate * info.theoretical_boost,
            theoretical_boost: info.theoretical_boost,
            shield_score: info.shield_primes.len(),
            shield_primes: join_primes(&info.shield_primes),
            trap_primes: join_primes(&info.trap_primes),
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

    fn approx_eq(a: f64, b: f64) -> bool {
        (a - b).abs() < 1e-9
    }

    #[test]
    fn local_effect_cases() {
        // Gap 2 is 2 = -1 (mod 3): twin primes force p = 2 (mod 3), so
        // S = 2p + 1 = 2 (mod 3). Shielded, not exposed.
        assert_eq!(local_effect(2, 3), LocalEffect::Shield);
        assert_eq!(local_effect(4, 3), LocalEffect::Shield);
        assert_eq!(local_effect(6, 3), LocalEffect::Trap);
        assert_eq!(local_effect(2, 5), LocalEffect::Neutral);
        assert_eq!(local_effect(4, 5), LocalEffect::Shield);
        assert_eq!(local_effect(34, 11), LocalEffect::Shield); // 34 = 3*11 + 1
        assert_eq!(local_effect(210, 7), LocalEffect::Trap);
    }

    #[test]
    fn local_factor_values() {
        assert!(approx_eq(local_factor(4, 3), 1.5));
        assert!(approx_eq(local_factor(6, 3), 0.75)); // the "Mod 3 Trap"
        assert!(approx_eq(local_factor(2, 5), 5.0 * 2.0 / (4.0 * 3.0)));
        assert!(approx_eq(local_factor(4, 5), 1.25));
        assert!(approx_eq(local_factor(10, 5), 5.0 * 3.0 / 16.0));
    }

    #[test]
    fn shield_and_trap_lists() {
        let info_2 = calculate_shielding_info(2);
        assert_eq!(info_2.shield_primes[..1], [3]);
        assert!(info_2.trap_primes.is_empty());

        let info_4 = calculate_shielding_info(4);
        assert_eq!(info_4.shield_primes[..2], [3, 5]);
        assert!(info_4.trap_primes.is_empty());

        let info_6 = calculate_shielding_info(6);
        assert_eq!(info_6.shield_primes[..2], [5, 7]);
        assert_eq!(info_6.trap_primes, [3]);

        // Gap 34: 33 * 35 = 3 * 11 * 5 * 7, so it is a quadruple shield.
        let info_34 = calculate_shielding_info(34);
        assert_eq!(info_34.shield_primes[..4], [3, 5, 7, 11]);
        assert_eq!(info_34.trap_primes, [17]); // 34 = 2 * 17, a mild trap

        // Gap 210 = 2 * 3 * 5 * 7 is trapped by three primes at once.
        let info_210 = calculate_shielding_info(210);
        assert_eq!(info_210.trap_primes, [3, 5, 7]);
        assert_eq!(info_210.shield_primes[..2], [11, 19]); // 209 = 11 * 19
    }

    #[test]
    fn gap_4_to_gap_2_ratio_is_exactly_three_halves() {
        // Gaps 2 and 4 are both shielded by 3 and both neutral for every
        // prime >= 7, so the only difference is the factor at 5:
        // 1.25 / (5*2/(4*3)) = 1.5. The 10^10 run observes 0.1439 / 0.0959.
        let r = calculate_shielding_info(4).theoretical_boost
            / calculate_shielding_info(2).theoretical_boost;
        assert!(approx_eq(r, 1.5));
    }

    #[test]
    fn boosts_are_ordered_as_observed() {
        let boost = |g| calculate_shielding_info(g).theoretical_boost;
        // Observed at 10^10: 34 (0.206) > 14 (0.168) > 4 (0.144) > 2 (0.096)
        // > 6 (0.090) > 30 (0.059) and 210 (0.065) below 2.
        assert!(boost(34) > boost(14));
        assert!(boost(14) > boost(4));
        assert!(boost(4) > boost(2));
        assert!(boost(2) > boost(6));
        assert!(boost(6) > boost(210));
        assert!(boost(210) > boost(30));
    }
}
