use bitvec::prelude::*;
use rayon::prelude::*;
use std::collections::VecDeque;
use std::time::{Duration, Instant};

#[derive(Clone, Copy, Debug)]
pub struct SieveTuning {
    pub min_chunk_words: usize,
    pub chunk_multiplier: usize,
}

#[derive(Clone, Copy, Debug)]
pub struct ChunkPlan {
    pub segment_words: usize,
    pub rayon_threads: usize,
    pub chunk_size_words: usize,
    pub approx_chunks_per_segment: usize,
}

#[derive(Clone, Copy, Debug)]
pub struct BenchmarkResult {
    pub segment_size_kb: usize,
    pub min_chunk_words: usize,
    pub chunk_multiplier: usize,
    pub rayon_threads: usize,
    pub chunk_size_words: usize,
    pub approx_chunks_per_segment: usize,
    pub measured_segments: usize,
    pub elapsed: Duration,
}

impl BenchmarkResult {
    pub fn segments_per_second(self) -> f64 {
        self.measured_segments as f64 / self.elapsed.as_secs_f64()
    }
}

fn effective_chunk_size_words(segment_words: usize, tuning: SieveTuning) -> usize {
    let num_threads = rayon::current_num_threads();
    let multiplier = tuning.chunk_multiplier.max(1);
    let min_chunk_words = tuning.min_chunk_words.max(1);

    if num_threads > 0 {
        segment_words
            .div_ceil(num_threads * multiplier)
            .max(min_chunk_words)
    } else {
        segment_words.max(1)
    }
}

pub fn describe_chunk_plan(segment_size_bytes: usize, tuning: SieveTuning) -> ChunkPlan {
    let segment_words = (segment_size_bytes / 8).max(1);
    let rayon_threads = rayon::current_num_threads();
    let chunk_size_words = effective_chunk_size_words(segment_words, tuning);
    let approx_chunks_per_segment = segment_words.div_ceil(chunk_size_words);

    ChunkPlan {
        segment_words,
        rayon_threads,
        chunk_size_words,
        approx_chunks_per_segment,
    }
}

/// An iterator that generates primes up to a given limit using a segmented sieve.
pub struct PrimeIterator {
    limit: u64,
    sqrt_limit: u64,
    base_primes: Vec<u32>,
    sieve_state: SieveState,
    segment_size_bits: u64,
    tuning: SieveTuning,
}

enum SieveState {
    Base(usize), // Index into base_primes
    Segmented {
        segment_start: u64,
        segment: BitVec<u64, Lsb0>,
        segment_index: usize,
    },
}

impl PrimeIterator {
    pub fn new(limit: u64, segment_size_bytes: usize, tuning: SieveTuning) -> Self {
        let sqrt_limit = (limit as f64).sqrt() as u64;

        let mut base_sieve = bitvec![u8, Lsb0; 1; (sqrt_limit + 1) as usize];
        base_sieve.set(0, false);
        base_sieve.set(1, false);

        for i in 2..=(sqrt_limit as f64).sqrt() as u64 {
            if base_sieve[i as usize] {
                for j in (i * i..=sqrt_limit).step_by(i as usize) {
                    base_sieve.set(j as usize, false);
                }
            }
        }

        let base_primes: Vec<u32> = base_sieve.iter_ones().map(|i| i as u32).collect();

        Self {
            limit,
            sqrt_limit,
            base_primes,
            sieve_state: SieveState::Base(0),
            segment_size_bits: (segment_size_bytes * 8) as u64,
            tuning,
        }
    }

    fn sieve_segment(
        start: u64,
        end: u64,
        base_primes: &[u32],
        tuning: SieveTuning,
    ) -> BitVec<u64, Lsb0> {
        let mut segment = bitvec![u64, Lsb0; 0; (end - start) as usize]; // 0 means prime

        // Parallelize over memory chunks (Domain Decomposition).
        // CRITICAL OPTIMIZATION: Thread-Aligned Chunking.
        // We calculate chunk_size such that we have roughly 1 chunk per thread.
        // This minimizes the overhead of calculating the initial prime offset (modulo),
        // which can be expensive if done for many small chunks.
        let raw_slice = segment.as_raw_mut_slice();
        let len_u64 = raw_slice.len();
        let chunk_size = effective_chunk_size_words(len_u64, tuning);

        // Optimization: Pre-mark even numbers.
        // If start is even: bit 0 is even (composite). Pattern 1010... (0x55...)
        // If start is odd: bit 0 is odd (prime?). Pattern 0101... (0xAA...)
        let pattern = if start.is_multiple_of(2) {
            0x5555555555555555
        } else {
            0xAAAAAAAAAAAAAAAA
        };

        raw_slice
            .par_chunks_mut(chunk_size)
            .enumerate()
            .for_each(|(chunk_idx, chunk)| {
                // Initialize chunk with even/odd pattern. Replaces loop for p=2.
                chunk.fill(pattern);

                // Determine the range of global bits this chunk covers
                let chunk_start_word_idx = chunk_idx * chunk_size;
                let chunk_start_bit = start + (chunk_start_word_idx as u64 * 64);
                let chunk_len_bits = chunk.len() as u64 * 64;

                for &p_u32 in base_primes {
                    if p_u32 == 2 {
                        continue;
                    } // Handled by pattern fill

                    let p = p_u32 as u64;
                    let p_sq = p * p;

                    // Calculate the bit offset within the chunk for the first multiple of p
                    let start_bit_in_chunk = if chunk_start_bit < p_sq {
                        // The first multiple we care about is p*p.
                        // Calculate offset of p*p relative to chunk_start_bit.
                        p_sq.saturating_sub(chunk_start_bit)
                    } else {
                        // chunk_start_bit >= p*p.
                        // Find the smallest k >= 0 such that (chunk_start_bit + k) is a multiple of p.
                        let rem = chunk_start_bit % p;
                        if rem == 0 {
                            0
                        } else {
                            p - rem
                        }
                    };

                    // Mark composites in this chunk
                    let mut local_bit_idx = start_bit_in_chunk;
                    while local_bit_idx < chunk_len_bits {
                        let word_idx = (local_bit_idx / 64) as usize;
                        let bit_idx = (local_bit_idx % 64) as usize;

                        chunk[word_idx] |= 1 << bit_idx;

                        local_bit_idx += p;
                    }
                }
            });

        if start == 0 {
            if !segment.is_empty() {
                segment.set(0, true); // 0 is not prime
            }
            if segment.len() > 1 {
                segment.set(1, true); // 1 is not prime
            }
            if segment.len() > 2 {
                segment.set(2, false); // 2 IS prime (was marked by pattern)
            }
        }

        segment
    }
}

impl Iterator for PrimeIterator {
    type Item = u64;

    fn next(&mut self) -> Option<Self::Item> {
        loop {
            match &mut self.sieve_state {
                SieveState::Base(index) => {
                    if *index < self.base_primes.len() {
                        let prime = self.base_primes[*index] as u64;
                        *index += 1;
                        if prime > self.limit {
                            return None;
                        }
                        return Some(prime);
                    } else {
                        let segment_start = self.sqrt_limit + 1;
                        let segment_end =
                            (segment_start + self.segment_size_bits).min(self.limit + 1);
                        let segment = Self::sieve_segment(
                            segment_start,
                            segment_end,
                            &self.base_primes,
                            self.tuning,
                        );
                        self.sieve_state = SieveState::Segmented {
                            segment_start,
                            segment,
                            segment_index: 0,
                        };
                    }
                }
                SieveState::Segmented {
                    segment_start,
                    segment,
                    segment_index,
                } => {
                    let raw = segment.as_raw_slice();
                    let start_idx = *segment_index;
                    let start_word = start_idx / 64;
                    let start_bit = start_idx % 64;

                    #[allow(clippy::needless_range_loop)]
                    for word_idx in start_word..raw.len() {
                        let mut word = raw[word_idx];

                        // If this is the first word, mask out bits we've already processed
                        if word_idx == start_word {
                            // Set bits 0..start_bit to 1 (treat as composite/processed)
                            // so trailing_ones() skips them.
                            if start_bit < 64 {
                                let mask = (1u64 << start_bit) - 1;
                                word |= mask;
                            } else {
                                word = u64::MAX;
                            }
                        }

                        // If word is not all ones, there is a zero (prime)
                        if word != u64::MAX {
                            let bit_idx = word.trailing_ones() as usize;
                            let found_index = word_idx * 64 + bit_idx;

                            *segment_index = found_index + 1;
                            let prime = *segment_start + found_index as u64;
                            if prime > self.limit {
                                return None;
                            }
                            return Some(prime);
                        }
                    }

                    // Segment exhausted
                    *segment_index = segment.len();

                    *segment_start += self.segment_size_bits;
                    if *segment_start > self.limit {
                        return None;
                    }
                    let segment_end = (*segment_start + self.segment_size_bits).min(self.limit + 1);
                    *segment = Self::sieve_segment(
                        *segment_start,
                        segment_end,
                        &self.base_primes,
                        self.tuning,
                    );
                    *segment_index = 0;
                }
            }
        }
    }
}

pub struct PrimalityChecker {
    limit: u64,
    sqrt_limit: u64,
    base_primes: Vec<u32>,
    known_primes_under_sqrt: BitVec<u8, Lsb0>,

    cached_segments: VecDeque<(u64, BitVec<u64, Lsb0>)>,
    cache_size: usize,
    segment_size_bits: u64,
    tuning: SieveTuning,
}

impl PrimalityChecker {
    pub fn new(limit: u64, segment_size_bytes: usize, tuning: SieveTuning) -> Self {
        let sqrt_limit = (limit as f64).sqrt() as u64;

        let mut base_sieve = bitvec![u8, Lsb0; 1; (sqrt_limit + 1) as usize];
        base_sieve.set(0, false);
        base_sieve.set(1, false);

        for i in 2..=(sqrt_limit as f64).sqrt() as u64 {
            if base_sieve[i as usize] {
                for j in (i * i..=sqrt_limit).step_by(i as usize) {
                    base_sieve.set(j as usize, false);
                }
            }
        }

        let base_primes: Vec<u32> = base_sieve.iter_ones().map(|i| i as u32).collect();

        Self {
            limit,
            sqrt_limit,
            base_primes,
            known_primes_under_sqrt: base_sieve,
            cached_segments: VecDeque::with_capacity(8),
            cache_size: 8,
            segment_size_bits: (segment_size_bytes * 8) as u64,
            tuning,
        }
    }

    pub fn ensure_range(&mut self, start: u64, end: u64) {
        if start > self.limit {
            return;
        }
        let effective_start = std::cmp::max(start, self.sqrt_limit + 1);
        if effective_start > end {
            return;
        }

        let start_seg = (effective_start / self.segment_size_bits) * self.segment_size_bits;
        let end_seg = (end / self.segment_size_bits) * self.segment_size_bits;

        let mut seg = start_seg;
        while seg <= end_seg {
            // Check if seg in cache
            if !self.cached_segments.iter().any(|(s, _)| *s == seg) {
                let segment_end = seg + self.segment_size_bits;
                let new_segment =
                    PrimeIterator::sieve_segment(seg, segment_end, &self.base_primes, self.tuning);

                if self.cached_segments.len() >= self.cache_size {
                    self.cached_segments.pop_front();
                }
                self.cached_segments.push_back((seg, new_segment));
            }
            seg += self.segment_size_bits;
        }
    }

    pub fn is_prime(&mut self, n: u64) -> bool {
        if n > self.limit {
            return false;
        }
        if n <= self.sqrt_limit {
            return self.known_primes_under_sqrt[n as usize];
        }

        let segment_start = (n / self.segment_size_bits) * self.segment_size_bits;

        for (start, segment) in &self.cached_segments {
            if *start == segment_start {
                let index = (n - segment_start) as usize;
                return !segment[index];
            }
        }

        let segment_end = segment_start + self.segment_size_bits;
        let new_segment = PrimeIterator::sieve_segment(
            segment_start,
            segment_end,
            &self.base_primes,
            self.tuning,
        );

        let is_p = !new_segment[(n - segment_start) as usize];

        if self.cached_segments.len() >= self.cache_size {
            self.cached_segments.pop_front();
        }
        self.cached_segments.push_back((segment_start, new_segment));

        is_p
    }
}

pub fn benchmark_sieve_configs(
    max_n: u64,
    segment_kb_list: &[usize],
    min_chunk_words_list: &[usize],
    chunk_multiplier: usize,
    warmup_segments: usize,
    measured_segments: usize,
) -> Vec<BenchmarkResult> {
    let mut results = Vec::new();
    let segment_kb_list: Vec<usize> = segment_kb_list.iter().copied().filter(|v| *v > 0).collect();
    let min_chunk_words_list: Vec<usize> = min_chunk_words_list
        .iter()
        .copied()
        .filter(|v| *v > 0)
        .collect();

    for &segment_size_kb in &segment_kb_list {
        let segment_size_bytes = segment_size_kb * 1024;
        for &min_chunk_words in &min_chunk_words_list {
            let tuning = SieveTuning {
                min_chunk_words,
                chunk_multiplier,
            };
            let iterator = PrimeIterator::new(max_n, segment_size_bytes, tuning);
            let mut start = iterator.sqrt_limit + 1;
            let mut checksum = 0u64;

            for _ in 0..warmup_segments {
                let end = start + iterator.segment_size_bits;
                let segment =
                    PrimeIterator::sieve_segment(start, end, &iterator.base_primes, tuning);
                checksum ^= segment.as_raw_slice().first().copied().unwrap_or(0);
                start = end;
            }

            let bench_start = Instant::now();
            for _ in 0..measured_segments {
                let end = start + iterator.segment_size_bits;
                let segment =
                    PrimeIterator::sieve_segment(start, end, &iterator.base_primes, tuning);
                checksum ^= segment.as_raw_slice().first().copied().unwrap_or(0);
                start = end;
            }
            let elapsed = bench_start.elapsed();

            let plan = describe_chunk_plan(segment_size_bytes, tuning);
            let _ = checksum;
            results.push(BenchmarkResult {
                segment_size_kb,
                min_chunk_words,
                chunk_multiplier,
                rayon_threads: plan.rayon_threads,
                chunk_size_words: plan.chunk_size_words,
                approx_chunks_per_segment: plan.approx_chunks_per_segment,
                measured_segments,
                elapsed,
            });
        }
    }

    results
}
