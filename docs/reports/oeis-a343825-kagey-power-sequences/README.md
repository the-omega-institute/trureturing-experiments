# OEIS A343825: least end of a sequence whose product is an n-th power

Check for the open-problem preregistration on Kagey's conjecture in OEIS A343825: T(p, k) = A071537(k) for fixed k and sufficiently large prime p.

Run: `python3 check.py 150` (Python 3, standard library only).

For a prime p, a sequence k = b_1 ≤ … ≤ b_t = m with multiplicities below p and a p-th-power product is the same as a multiplicity vector c_j ∈ [0, p − 1] (j = k, …, m). It needs c_k ≥ 1, c_m ≥ 1 and Σ c_j v(j) ≡ 0 (mod p), where v(j) is the exponent vector of j. The program decides this by Gaussian elimination over F_p, and decides A071537 by elimination over Q.

It checks:
- the published data of A071537;
- the prime rows 2, 3, 5 and 7 of the A343825 table;
- T(p, k) = A071537(k) for 1 ≤ k ≤ 150 at p = 31, 37, 41 and 43. At p = 37 the exceptions k = 144, 147, 150 occur; they agree again at 41 and 43, as "sufficiently large" allows.

`run.log` holds the output of the run cited by the preregistration.
