use bitvec::prelude::*;
use rayon::prelude::*;
use std::collections::VecDeque;
use std::sync::RwLock;

/// An iterator that generates primes up to a given limit using a segmented sieve.
pub struct PrimeIterator {
    limit: u64,
    sqrt_limit: u64,
    base_primes: Vec<u32>,
    sieve_state: SieveState,
    segment_size_bits: u64,
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
    pub fn new(limit: u64, segment_size_bytes: usize) -> Self {
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
        }
    }

    fn sieve_segment(start: u64, end: u64, base_primes: &[u32]) -> BitVec<u64, Lsb0> {
        let mut segment = bitvec![u64, Lsb0; 0; (end - start) as usize]; // 0 means prime

        // Parallelize over memory chunks (Domain Decomposition).
        // CRITICAL OPTIMIZATION: Thread-Aligned Chunking.
        // We calculate chunk_size such that we have roughly 1 chunk per thread.
        // This minimizes the overhead of calculating the initial prime offset (modulo),
        // which can be expensive if done for many small chunks.
        let raw_slice = segment.as_raw_mut_slice();
        let len_u64 = raw_slice.len();

        let num_threads = rayon::current_num_threads();
        let multiplier = 4; // Create 4x chunks to allow P-cores to steal work from E-cores
        let chunk_size = if num_threads > 0 {
            len_u64.div_ceil(num_threads * multiplier).max(4096)
        } else {
            len_u64.max(1)
        };

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
                        let segment =
                            Self::sieve_segment(segment_start, segment_end, &self.base_primes);
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
                    *segment = Self::sieve_segment(*segment_start, segment_end, &self.base_primes);
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

    cached_segments: RwLock<VecDeque<(u64, BitVec<u64, Lsb0>)>>,
    cache_size: usize,
    segment_size_bits: u64,
}

impl PrimalityChecker {
    pub fn new(limit: u64, segment_size_bytes: usize) -> Self {
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
            cached_segments: RwLock::new(VecDeque::with_capacity(4)),
            cache_size: 4,
            segment_size_bits: (segment_size_bytes * 8) as u64,
        }
    }

    pub fn is_prime(&self, n: u64) -> bool {
        if n > self.limit {
            return false;
        }
        if n <= self.sqrt_limit {
            return self.known_primes_under_sqrt[n as usize];
        }

        let segment_start = (n / self.segment_size_bits) * self.segment_size_bits;

        // 1. Try to find in cache with Read Lock
        {
            let cache_read = self.cached_segments.read().unwrap();
            for (start, segment) in cache_read.iter() {
                if *start == segment_start {
                    let index = (n - segment_start) as usize;
                    return !segment[index];
                }
            }
        }

        // 2. Not found, acquire Write Lock
        let mut cache_write = self.cached_segments.write().unwrap();

        // 3. Double-check (another thread might have added it)
        for (start, segment) in cache_write.iter() {
            if *start == segment_start {
                let index = (n - segment_start) as usize;
                return !segment[index];
            }
        }

        // 4. Generate new segment
        let segment_end = segment_start + self.segment_size_bits;
        let new_segment =
            PrimeIterator::sieve_segment(segment_start, segment_end, &self.base_primes);

        let is_p = !new_segment[(n - segment_start) as usize];

        if cache_write.len() >= self.cache_size {
            cache_write.pop_front();
        }
        cache_write.push_back((segment_start, new_segment));

        is_p
    }
}
