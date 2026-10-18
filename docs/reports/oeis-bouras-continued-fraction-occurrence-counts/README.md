# OEIS A372761, A372763, A363347, A363102: values and multiplicities of Bouras' continued-fraction denominators

The sequences are the denominators of K_n(x) = 1/(2 − 3/(3 − 4/(4 − ⋯ (n−1) − n/x))) for n ≥ 3, with x = n+4 (A372761), x = n+5 (A372763), x = −4 (A363347) and x = −2 (A363102). The statements checked are:
- **A372761, Conjecture 1.** Except for 4, every denominator is 1 or a prime.
- **A372763, Conjecture 2.** Except for 2 and 3, every prime occurs exactly once. It is refuted: 11 occurs at n = 10 and n = 21.
- **A363347, Conjecture 3.** Every prime value other than 5 occurs exactly twice.
- **A363102, Conjecture 2.** Every prime value occurs exactly twice or three times.

`check.py` (Python 3 with sympy) does three things:
1. **Exact evaluation.** It evaluates the nested fraction exactly for 3 ≤ n ≤ N1, compares with the published data, and checks the closed form and the denominator formula den = |D|/gcd(|D|, n!/2), D = (n−1)x − n(n−2).
2. **Valuation formula.** It computes the denominators for n ≤ N2 from the p-adic deficits of |D| against n!/2.
3. **Statements and predictions.** It checks the four statements and the proofs' predictions for every prime p ≤ N2/2. All occurrences of such a prime lie below 2p.

Run `python3 check.py 300 20000`. Exit 0 and `RESULT PASS` mean every check holds.
- **A363102.** The primes occurring three times up to 10000 are 7, 41 and 239, at the third index n = √(2p²+2) = 10, 58, 338.
- **A372763.** The only prime occurring twice is 11.
