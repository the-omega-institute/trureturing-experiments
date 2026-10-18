# OEIS A363558, A363559, A363561, A363569: Hanna's bilateral sums Σ x^n (c + x^n)^(2n)

For c = 2, 3, √2, i, P. D. Hanna conjectures:
- a(2^n) = c^(2^(n+1)) for n > 0;
- a(p) = p·c^(p−1) + c^(2p) for primes p > 3.

**Coefficient formula.** The term n ≥ 1 contributes C(2n, m/n − 1)·c^(2n−m/n+1) to x^m when n | m. The term n = −k contributes C(−2k, m/k − 2k + 1)·c^(m/k−2k+1) when k | m.

**The prime case.** For a prime p only two terms remain.

**The case m = 2^N.** The positive term at n = 2^i (i < N) cancels the negative term at k = 2^(N−1−i): same power of c, opposite coefficient, because C(−a, b) = (−1)^b C(a+b−1, b) with b odd. This leaves c^(2^(N+1)).

`check.py` (Python 3 with sympy, c symbolic) checks:
- the formula against a direct expansion of the bilateral sum;
- the published data;
- both statements, as polynomials in c and for the four values of c.

Run `python3 check.py 300 12`. Exit 0 and `RESULT PASS` mean every check holds.
