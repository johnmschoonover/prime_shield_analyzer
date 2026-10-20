# OEIS A401296

`a(n)` is the smallest even `k > 0` with `k == +1` or `-1 (mod p)` for each of
the first `n` odd primes; equivalently the smallest even square root of unity
modulo the product of the first `n` odd primes. These are the smallest prime
gaps "shielded" against every odd prime up to `prime(n+1)` in the sense of the
analyzer in this repository.

- `a401296.py` computes the sequence by meet-in-the-middle over the `2^n` sign
  patterns, `O(2^(n/2) log)` instead of `O(2^n)`. `python3 a401296.py 40`
  prints `n a(n)` for `n = 1..40` in about twenty seconds.
- `b401296.txt` is the resulting b-file, `n = 1..40`, verified three ways:
  against brute force for `n <= 18`, against an independent computation by an
  OEIS editor for `n = 19..21`, and against a rerun with a different split
  point for `n <= 30`. Every term is even, satisfies `k^2 == 1 (mod P_n)`, and
  is at most `P_n - 1`.
