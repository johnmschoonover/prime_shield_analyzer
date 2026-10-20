# `prime-shield-analyzer`

A Rust tool for measuring a structural bias in sums of consecutive primes. For each pair of consecutive primes $(p_n, p_{n+1})$ it asks whether

$$S_n = p_n + p_{n+1} - 1$$

is itself prime, and tallies the hit rate by the prime gap $g = p_{n+1} - p_n$. The hit rate is far from uniform: at $N = 10^{12}$ it ranges from 5.6% (gap 210) to 17.1% (gap 34). This repository contains the sieve that measures it, the model that predicts it, and an interactive report comparing the two.

> **[View the interactive report](https://johnmschoonover.github.io/prime_shield_analyzer/)**

The write-ups live on [theschoonover.net](https://theschoonover.net/writing): *Defeating the Sieve* (Parts 1 and 2), and Part 3, which corrects the model used in the first two.

---

## The model: shields, traps, and neutral primes

Both parents are prime, so for any odd prime $q$ smaller than them, $p_n$ avoids two residues: $p_n \not\equiv 0$ and $p_n \not\equiv -g \pmod q$. The sum $S = 2p_n + g - 1$ is divisible by $q$ exactly when $p_n \equiv (1-g)/2 \pmod q$. Comparing that bad residue with the two the parents already exclude gives three cases.

| Case | Condition | $S$ divisible by $q$ for | Factor vs. a random odd integer |
|---|---|---|---|
| **Shield** | $g \equiv \pm 1 \pmod q$ | no residue | $\dfrac{q}{q-1}$ |
| **Trap** | $g \equiv 0 \pmod q$ | 1 of $q-1$ residues | $\dfrac{q(q-2)}{(q-1)^2}$ |
| **Neutral** | otherwise | 1 of $q-2$ residues | $\dfrac{q(q-3)}{(q-1)(q-2)}$ |

The predicted hit rate for a gap is the size-matched baseline, the mean of $2/\ln S$ over its pairs, multiplied by the product of these factors over all odd primes (the tool truncates at 1000; the tail is negligible). This is the Bateman–Horn local factor for the linear system $(n,\; n+g,\; 2n+g-1)$, conditioned on the first two forms being prime, so none of it is new mathematics. What the tool adds is the measurement.

Some consequences:

- **Gap 2 is shielded against 3.** Twin primes above 3 are always $p \equiv 2 \pmod 3$, so $S = 2p+1 \equiv 2 \pmod 3$. The first two write-ups and earlier versions of this README got this wrong.
- **Gap 4 beats gap 2 by exactly 3/2.** Both are shielded at 3 and neutral at every prime from 7 up. The only difference is at 5, where gap 4 is shielded and gap 2 is neutral: $1.25 / (10/12) = 1.5$. Observed at $10^{12}$: $0.1194 / 0.0796 = 1.4999$; at $10^{10}$: $0.1439 / 0.0959 = 1.5003$.
- **Gap 34 is a quadruple shield**, not a triple: $33 \cdot 35 = 3 \cdot 5 \cdot 7 \cdot 11$. In general a gap is shielded by every odd prime up to $p_k$ exactly when $g^2 - 1$ is divisible by $3 \cdot 5 \cdots p_k$.
- **The trap at 3 is a 0.75 penalty**, which is why gaps 6, 12, 30 and 210 sit below gap 2 despite having shields elsewhere.
- **Fit.** Over the 121 gaps with at least 100,000 pairs at $10^{12}$, predicted vs. observed rate has $R^2 = 0.997$ (0.996 over 59 gaps at $10^{10}$). The earlier one-sided rule (shield only at $g \equiv 1 \pmod 3$ and $g \equiv -1 \pmod q$) scores $R^2 \approx 0.30$ on the same data.
- **Residuals shrink with $N$.** The largest relative misses (gap 12 at +6.7%, gap 160 at +5.4% at $10^{12}$) fall monotonically across the $10^8$, $10^{10}$ and $10^{12}$ runs, as expected for Hardy–Littlewood secondary terms; the mean absolute relative miss over well-sampled gaps goes 3.1% → 2.6% → 2.3%.

This is not a way to find primes. Everything the gap tells you about $S$ is already visible in $S \bmod 3, 5, 7, \ldots$, which trial division reads off directly and to a far greater depth.

## Technology

- Rust, with a parallel segmented Sieve of Eratosthenes over `bitvec` (via `rayon`). $N = 10^{10}$ runs in about eight minutes on four cloud cores; $N = 10^{12}$ (37.6 billion primes) took about 3.3 hours on a desktop with `--segment-size-kb 512`. Memory stays under 100 MB.
- Per-gap tallies of pairs, successes, and the $2/\ln S$ baseline, written as CSV.
- A self-contained HTML report (Chart.js) with a predicted-vs-observed scatter, an oscillation series across the number line, and the gap spectrum.

## Installation and usage

```bash
git clone https://github.com/johnmschoonover/prime_shield_analyzer.git
cd prime_shield_analyzer
cargo build --release
```

Contributors can install the pre-commit hooks, which run `cargo fmt`, `cargo clippy`, `cargo check`, and `cargo test`:

```bash
pip install pre-commit
pre-commit install
```

Run the analyzer with `--max-exponent E` for $N = 10^E$:

```bash
./target/release/prime_shield_analyzer \
    --max-exponent 8 \
    --gaps 2,4,6,14,30,34,210 \
    --web-report
```

| Argument | Meaning |
|---|---|
| `--max-exponent <E>` | Analyze primes up to $10^E$. |
| `--gaps <GAPS>` | Comma-separated even gaps to track in the oscillation series. Default `2,4,6,12,30`. |
| `--bins <N>` | Number of bins for the oscillation series. Default 1000. |
| `--segment-size-kb <KB>` | Sieve segment size, for cache tuning. Default 128. |
| `--output-dir <DIR>` | Output directory. Default `results`. |
| `--web-report` | Also write `index.html` and the JSON files it reads. |

## Output

`gap_spectrum.csv` has one row per gap size:

| Column | Meaning |
|---|---|
| `count`, `successes`, `success_rate` | Pairs seen, pairs whose $S$ was prime, and their ratio. |
| `baseline_rate` | Mean of $2/\ln S$ over the pairs: a random odd integer's chance. |
| `theoretical_boost` | Product of the local factors above. |
| `predicted_rate` | `baseline_rate × theoretical_boost`. |
| `shield_score`, `shield_primes`, `trap_primes` | Number and list of shielding primes, and the list of trapping primes. |

`oscillation_series.csv` gives per-bin counts and per-gap rates across the number line. `global_stats.csv` gives the overall totals. The HTML report's first chart plots `predicted_rate` against `success_rate` for every gap with at least 1,000 samples, with a dashed $y = x$ line and the $R^2$ in the header. The shipped `report/` directory is the $10^{12}$ run.

`scripts/summarize_spectrum.py` turns any `gap_spectrum.csv` into the $R^2$ figures (new model and old rule, at two sample thresholds), the largest residual, the gap 4 / gap 2 ratio, and a Markdown table for a chosen list of gaps:

```bash
python3 scripts/summarize_spectrum.py report/gap_spectrum.csv 2,4,6,14,30,34,210
```

---

## Python tool: `prime_shields.py`

`prime_shields.py` computes a sequence of gap sizes that are fully shielded under the one-sided rule: the smallest even $G_k$ with $G_k \equiv 1 \pmod 3$ and $G_k \equiv -1 \pmod{p}$ for every prime $5 \le p \le p_k$. It uses a ratchet: each step advances by the product of the primes already satisfied, which preserves their residues, until the next prime's residue is hit. The first terms are 4, 4, 34, 1924, 25024, 85084.

Under the symmetric $\pm 1$ rule the fully shielded gaps are smaller (2, 4, 34, 34, 274, 9164, 83656, 937364), and the smallest solution is not always an extension of the previous one, so that sequence needs a search over sign patterns rather than a ratchet. The symmetric sequence is OEIS [A401296](https://oeis.org/A401296), submitted October 2026; it is the even counterpart of [A323746](https://oeis.org/A323746). The one-sided variant is not in the OEIS.

```bash
pip install sympy
python3 prime_shields.py 6
```
