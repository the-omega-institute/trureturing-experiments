# Benfield–Lippard Question 4.34: exact corroboration

Source: Brennan Benfield and Oliver Lippard, *Connecting Zeros in Pisano Periods
to Prime Factors of K-Fibonacci Numbers*, arXiv:2407.20048v2, Question 4.34
and the period table following Theorem 4.28.

Run with Python 3 and no third-party packages:

```sh
python3 docs/reports/benfield-lippard-q434/check.py
```

The sequence is U₀ = 0, U₁ = 1, Uₙ₊₂ = K Uₙ₊₁ + Uₙ. The program independently
counts zeros in the shortest pair cycle, finds the least positive zero index,
and verifies their product is the period. A separate binary matrix-power
implementation checks the period and its minimality by its prime divisors.
All evidence uses integers; no floating-point calculation is used.

For positive even K and m = 2ᵃs with a ≥ 2 and odd s, the proposed answer sets
r = v₂(K), t = max(1, a − r + 1), c = v₂(π_K(s)), with π_K(1) = 1, and asserts
ω_K(m) = 2 if c > t and 1 otherwise. When ω_K(m) = 1, it further asserts
4 does not divide π_K(m) exactly when a ≤ r.

The saved output in `check.out` reports successful checks of:

- All 16 entries of the source's dyadic period table and three stated examples.
- The dyadic formula for 640 cases: even K from 2 to 128 and a from 1 to 10.
- Both clauses of the proposed answer on 12,288 pairs: even K from 2 to 96,
  every positive m ≤ 1024 divisible by 4. Orders 1 and 2 occur 6,607 and
  5,681 times respectively.
- Both clauses on 11,200 structured triples: even K from 2 to 64, a from 2 to 8,
  and odd s from 1 to 99. Orders 1 and 2 occur 8,791 and 2,409 times respectively.
- The rank and period CRT formulas and odd-component alternatives on both grids.
- The integer valuation formula v₂(Uₙ) = v₂(K) + v₂(n) − 1 for even n, and
  oddness for odd n, on 8,192 pairs: even K from 2 to 128 and n from 1 to 128.

The two main grids overlap and are reported separately. Their results are finite
corroboration of a proposed general argument, not a proof of the universal claim.
The program raises an exception on failure and prints `ALL_CHECKS_OK` on success.
