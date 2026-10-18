# OEIS A350532: perfect powers over GF(2) by number of terms

A350532's T(n,k) counts the degree-n polynomials over GF(2) of the form f(x)^m with m > 1 and exactly k nonzero terms. Peter Kagey's entry makes three conjectures:
- (1) T(n, n+1) = 1 iff n is a Mersenne prime;
- (2) T(2n, 2) = n;
- (3) T(2n, 3) = (n² − n)/2.

`check.py` (Python 3 with sympy) tests "f^m for some m > 1" as "the multiplicities in the GF(2) factorization have gcd > 1". It checks:
- the published rows n ≤ 10, by enumerating all polynomials;
- T(n, n+1) for n ≤ 70, where only the all-ones polynomial has n + 1 terms;
- T(2n, 2) for n ≤ 40, and T(n, 3) for n < 40.

Result:
- Conjecture (1) fails. T(n, n+1) = 1 holds exactly for n = 3, 7, 15, 31, 63, and the counterexample is 1 + x + ⋯ + x^15 = (1 + x)^15 = ((1 + x)^5)^3, with 15 not prime.
- Conjectures (2) and (3) hold in range.

Run it with `python3 check.py`. Exit 0 and `RESULT PASS` mean every check holds. `run.log` is the output (6 s).
