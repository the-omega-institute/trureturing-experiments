# OEIS A369378: Ordowski's conjecture that 2^(2^n+1) − 1 + 2^k is composite for k < 2^n

A369378(n) is the smallest k such that 2^(2^n+1) − 1 + 2^k is prime. T. Ordowski conjectures that for n > 0 the number is not prime for any k < 2^n. A Fermat-number covering proves it:
- **k = 0.** N = 2^(2^n+1) is not prime.
- **0 < k < 2^n.** Let j = v₂(k) < n. Then F_j = 2^(2^j) + 1 divides 2^(2^n+1) − 1 + 2^k, and N > F_j.

`check.py` (Python 3 with sympy) checks the covering for every k < 2^n, 1 ≤ n ≤ NMAX. It also recomputes the published terms (with `--full` including a(4) = 1968). Finally it checks that k = 2^n gives a prime for n = 1, 2, so the bound k < 2^n is sharp there.

Run `python3 check.py 12 --full`. Exit 0 and `RESULT PASS` mean every check holds.
