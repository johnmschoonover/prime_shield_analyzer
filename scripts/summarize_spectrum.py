#!/usr/bin/env python3
"""Summarize a gap_spectrum.csv from prime_shield_analyzer.

Prints R^2 of predicted vs observed success rate at two sample thresholds,
the R^2 the old one-sided rule would have scored on the same data, and a
Markdown table for a chosen list of gaps. Standard library only.

Usage: summarize_spectrum.py path/to/gap_spectrum.csv [gap,gap,...]
"""
import csv
import sys

OLD_PRIMES = [5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61,
              67, 71, 73, 79, 83, 89, 97]


def old_rule_boost(g):
    b = 1.5 if g % 3 == 1 else 1.0
    for q in OLD_PRIMES:
        if g % q == q - 1:
            b *= q / (q - 1)
    return b


def r_squared(xs, ys):
    n = len(xs)
    mx, my = sum(xs) / n, sum(ys) / n
    sxy = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
    sxx = sum((x - mx) ** 2 for x in xs)
    syy = sum((y - my) ** 2 for y in ys)
    return sxy * sxy / (sxx * syy)


def main(path, gaps):
    with open(path, newline="") as f:
        rows = [r for r in csv.DictReader(f) if int(r["gap_size"]) % 2 == 0]
    by_gap = {int(r["gap_size"]): r for r in rows}

    for thr in (1_000, 100_000):
        big = [r for r in rows if int(r["count"]) >= thr]
        obs = [float(r["success_rate"]) for r in big]
        pred = [float(r["predicted_rate"]) for r in big]
        old = [old_rule_boost(int(r["gap_size"])) for r in big]
        print(f"gaps with >= {thr:>7,} pairs: {len(big):>3}   "
              f"R^2 new = {r_squared(pred, obs):.4f}   "
              f"R^2 old rule = {r_squared(old, obs):.4f}")

    big = [r for r in rows if int(r["count"]) >= 100_000]
    worst = max(big, key=lambda r: abs(float(r["success_rate"]) - float(r["predicted_rate"])))
    print(f"largest miss (>= 100k pairs): gap {worst['gap_size']} "
          f"observed {float(worst['success_rate']):.4f} predicted {float(worst['predicted_rate']):.4f}")
    g2, g4 = by_gap[2], by_gap[4]
    print(f"gap 4 / gap 2: observed {float(g4['success_rate']) / float(g2['success_rate']):.4f}, "
          f"predicted {float(g4['predicted_rate']) / float(g2['predicted_rate']):.4f}")
    print()

    print("| Gap | Pairs | Observed | Predicted | Shields | Traps |")
    print("|---|---|---|---|---|---|")
    for g in gaps:
        r = by_gap.get(g)
        if r is None:
            print(f"| {g} | 0 | – | – | – | – |")
            continue
        shields = ", ".join(r["shield_primes"].split(",")[:4]) or "none"
        traps = ", ".join(r["trap_primes"].split(",")) or "none"
        print(f"| {g} | {int(r['count']):,} | {float(r['success_rate']):.4f} | "
              f"{float(r['predicted_rate']):.4f} | {shields} | {traps} |")


if __name__ == "__main__":
    if len(sys.argv) < 2:
        sys.exit(__doc__)
    gap_list = [int(x) for x in sys.argv[2].split(",")] if len(sys.argv) > 2 else [2, 4, 6, 8, 12, 14, 30, 34, 210]
    main(sys.argv[1], gap_list)
