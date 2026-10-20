#!/usr/bin/env python3
"""A401296 via meet-in-the-middle.

a(n) = smallest even k > 0 with k == +-1 (mod p) for each of the first n odd
primes. Split the primes into A (with the modulus 2, residue 0, folded in) and
B. For each sign pattern on A compute k_A mod 2*P_A; for each pattern on B
compute r_B mod P_B. Any k with k == k_A (mod 2P_A) is k_A + 2P_A*t, and the
B-condition fixes t == (r_B - k_A) * inv(2P_A) (mod P_B). The smallest such k
over all r_B is a successor query in a sorted table of r_B*inv mod P_B, so the
whole search costs O(2^(n/2) log) instead of O(2^n).
"""
import bisect
import sys
from itertools import product


def odd_primes(n):
    ps, c = [], 3
    while len(ps) < n:
        if all(c % p for p in ps):
            ps.append(c)
        c += 2
    return ps


def residues(primes, moduli_extra=()):
    """All CRT residues of +-1 patterns over `primes` (plus fixed extras)."""
    mods = list(moduli_extra) + [(p, (1, p - 1)) for p in primes]
    M = 1
    for m, _ in mods:
        M *= m
    out = [0]
    cur = 1
    for m, choices in mods:
        inv = pow(cur, -1, m)
        nxt = []
        for r in out:
            for c in choices:
                nxt.append(r + cur * (((c - r) * inv) % m))
        out = nxt
        cur *= m
    return out, M


def a(n):
    ps = odd_primes(n)
    half = n // 2
    A, B = ps[:half], ps[half:]
    ra, MA = residues(A, moduli_extra=[(2, (0,))])   # MA = 2 * P_A, all even
    rb, MB = residues(B)
    inv = pow(MA, -1, MB)
    table = sorted((r * inv) % MB for r in rb)
    best = None
    for k_a in ra:
        c = (k_a * inv) % MB
        i = bisect.bisect_left(table, c)
        s = table[i] if i < len(table) else table[0] + MB
        k = k_a + MA * ((s - c) % MB)
        if k == 0:
            k = MA * MB
        if best is None or k < best:
            best = k
    return best


if __name__ == "__main__":
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 21
    for n in range(1, N + 1):
        print(n, a(n), flush=True)
