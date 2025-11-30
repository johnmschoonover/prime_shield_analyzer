use serde::Serialize;
use std::collections::HashMap;

#[derive(Debug, Serialize)]
pub struct BinStats {
    pub bin_start: u64,
    pub bin_end: u64,
    pub prime_count_p: u64,
    pub prime_count_s: u64,
    #[serde(skip)]
    pub gap_successes: HashMap<u64, u64>, // Map<GapSize, Successes>
    #[serde(skip)]
    pub gap_occurrences: HashMap<u64, u64>, // Map<GapSize, Occurrences>
}

impl BinStats {
    fn new(start: u64, end: u64, target_gaps: &[u64]) -> Self {
        let mut gap_successes = HashMap::new();
        let mut gap_occurrences = HashMap::new();
        for &g in target_gaps {
            gap_successes.insert(g, 0);
            gap_occurrences.insert(g, 0);
        }
        Self {
            bin_start: start,
            bin_end: end,
            prime_count_p: 0,
            prime_count_s: 0,
            gap_successes,
            gap_occurrences,
        }
    }
}

#[derive(Debug)]
pub struct Statistics {
    pub total_primes: u64,
    pub total_s_primes: u64,
    pub gap_spectrum: HashMap<u64, (u64, u64)>, // Map<GapSize, (Occurrences, Successes)>
    pub bins: Vec<BinStats>,
    bin_size: u64,
    max_n_analysis_range: u64,
    pub target_gaps: Vec<u64>, // Store this for output.rs
}

impl Statistics {
    pub fn new(max_n: u64, num_bins: usize, target_gaps: &[u64]) -> Self {
        let max_n_analysis_range = max_n * 2;
        let bin_size = (max_n_analysis_range as f64 / num_bins as f64).ceil() as u64;

        let bins = (0..num_bins)
            .map(|i| {
                let start = (i as u64) * bin_size;
                let end = start + bin_size - 1;
                BinStats::new(start, end.min(max_n_analysis_range), target_gaps)
            })
            .collect();

        Self {
            total_primes: 0,
            total_s_primes: 0,
            gap_spectrum: HashMap::new(),
            bins,
            bin_size,
            max_n_analysis_range,
            target_gaps: target_gaps.to_vec(),
        }
    }

    pub fn get_bin_index(&self, n: u64) -> Option<usize> {
        if n > self.max_n_analysis_range {
            return None;
        }
        let index = (n / self.bin_size) as usize;
        if index < self.bins.len() {
            Some(index)
        } else {
            Some(self.bins.len() - 1)
        }
    }
}
