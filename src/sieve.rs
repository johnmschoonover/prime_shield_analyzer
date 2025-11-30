use bitvec::prelude::*;
use rayon::prelude::*;
use std::collections::VecDeque;
use std::sync::atomic::{AtomicU64, Ordering};

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

        // Unsafe cast to atomic slice. This is safe because u64 and AtomicU64 have the
        // same memory representation, and we are only performing atomic operations.
        let atomic_segment: &[AtomicU64] = unsafe {
            std::slice::from_raw_parts(
                segment.as_raw_slice().as_ptr() as *const AtomicU64,
                segment.as_raw_slice().len(),
            )
        };

        // Parallelize the marking of composites
        base_primes.par_iter().for_each(|&p_u32| {
            let p = p_u32 as u64;
            let mut mark_start = start.div_ceil(p) * p;
            if p * p > start {
                mark_start = p * p;
            }

            for i in (mark_start..end).step_by(p as usize) {
                let idx = (i - start) as usize;
                let word_idx = idx / 64;
                let bit_in_word = idx % 64;
                if word_idx < atomic_segment.len() {
                    atomic_segment[word_idx].fetch_or(1 << bit_in_word, Ordering::Relaxed);
                }
            }
        });

        if start == 0 {
            if !segment.is_empty() {
                segment.set(0, true);
            }
            if segment.len() > 1 {
                segment.set(1, true);
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
                    while *segment_index < segment.len() {
                        if !segment[*segment_index] {
                            let prime = *segment_start + *segment_index as u64;
                            *segment_index += 1;
                            if prime > self.limit {
                                return None;
                            }
                            return Some(prime);
                        }
                        *segment_index += 1;
                    }

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

    cached_segments: VecDeque<(u64, BitVec<u64, Lsb0>)>,
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
            cached_segments: VecDeque::with_capacity(4),
            cache_size: 4,
            segment_size_bits: (segment_size_bytes * 8) as u64,
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
        let new_segment =
            PrimeIterator::sieve_segment(segment_start, segment_end, &self.base_primes);

        let is_p = !new_segment[(n - segment_start) as usize];

        if self.cached_segments.len() >= self.cache_size {
            self.cached_segments.pop_front();
        }
        self.cached_segments.push_back((segment_start, new_segment));

        is_p
    }
}
