# De Castro's binomial formula for the p-adic order

Checks for the open-problem preregistrations on OEIS A343249–A343253 and A345531 and on Conjectures 1 and 2 of D. T. de Castro, "p-Adic order of positive integers via binomial coefficients", INTEGERS 22 (2022), #A61.

`L_p(k0, n)` is the lowest common denominator of `binomial(n, k)/n` over `0 < k <= k0` with `p ∤ k`. It divides `n`, so the formula `v_p(n) = log_p(n / L)` holds iff `n / L = p^v_p(n)`. The program uses exact integer arithmetic only.

Run: `python3 check.py 400 200` (Python 3, standard library only).

It checks:
- **A343249–A343253.** The published data, and for n ≤ 400 that the least k0 is the greatest power of a prime different from p dividing n (p = 2, 3, 5, 7, 11).
- **The valuation lemma.** For n ≤ 400 and every k, the r-part of the denominator of `binomial(n, k)/n` is at most `r^min(v_r(n), v_r(k))`, and it equals `r^v_r(n)` at `k = r^v_r(n)`.
- **INTEGERS Conjecture 1.** `L(S^min) = L(S^max) = n / p^v_p(n)` for primes p ≤ 13 and n ≤ 400. Every admissible set lies between these two, so the formula holds for all of them in that range.
- **INTEGERS Conjecture 2.** For composite q ≤ 30, primes `p_j | q` and c ≤ 2: `L(S^max) = 1` at `n = q^c p_j`, `v_q(n) = c`, and `n ≠ q^c`.
- **A345531.**
  - The published data.
  - For every composite q ≤ 200, the first n at which no k0 ∈ [1, n] works is `q * a(i)`, with `a(i)` the smallest prime power above the smallest prime factor of q.
  - The reading "k0 = n only" fails first at n = q.

The search for k0 stops as soon as L does not divide the target. This is exact, because `L(k0)` divides `L(k0 + 1)`.

`run.log` holds the output of the run cited by the preregistrations.
