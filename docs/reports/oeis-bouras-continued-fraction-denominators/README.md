# OEIS A356360, A370726, A372763, A362086, A363482: denominators of Bouras' continued fractions

Mohammed Bouras' five sequences are the denominators of K_n(x) = 1/(2 − 3/(3 − 4/(4 − ⋯ (n−1) − n/x))) for n ≥ 3:
- x = n+1, n+3, n+5 (A356360, A370726, A372763). The conjecture is that the denominators are 1 or primes.
- x = −3 and −5 (A362086, A363482). The conjecture is the same except for the single value 9, respectively 49.

`check.py` (Python 3 with sympy) does three things:
1. It evaluates the nested fraction exactly for 3 ≤ n ≤ N1, checking that no intermediate quotient is zero, and compares the result with the published data.
2. It checks the proof's closed form K_n(x) = (a(n−3)x − n·a(n−4))/((n−1)x − n(n−2)), where a = A051403, and the denominator formula den = |D|/gcd(|D|, n!/2) with D = (n−1)x − n(n−2).
3. It checks the conjectures with that formula for n ≤ N2, listing every composite denominator.

Run it with `python3 check.py 300 4000`. Exit 0 and `RESULT PASS` mean every check holds. The only composite denominators up to 4000 are 9 (A362086, n = 5) and 49 (A363482, n = 6).
