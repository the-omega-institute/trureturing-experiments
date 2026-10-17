# OEIS A357127: the two conjectures on the greatest prime factor of n^3 − 1

The sequence is a(n) = A081257(n) if A081257(n) > n, otherwise a(n) = 1 (n ≥ 2), where A081257(n) is the greatest prime factor of n³ − 1. The entry's two conjectures are:
- **Conjecture 1.** If a(n) ≠ 1, then a(n) = a(a(n) − n − 1).
- **Conjecture 2.** If n ≠ m and a(n) = a(m), then a(n) = gcd(n²+n+1, m²+m+1) = n + m + 1. It is refuted: a(35) = a(61) = 97, but gcd(1261, 3783) = 1261 = 13·97.

`check.py` (Python 3 with sympy) does four things for 2 ≤ n ≤ N:
- it computes a(n) by factorization;
- it checks Conjecture 1, computing the partner index directly;
- it lists the pairs failing Conjecture 2 and checks the predicted gcd formula gcd(n²+n+1, m²+m+1) = p·gcd((n²+n+1)/p, m − n);
- it checks the entry's comment that every prime value occurs exactly twice.

Run `python3 check.py 20000`. Exit 0 and `RESULT PASS` mean every check holds. The first failing pair, ordered by n, is (35, 61), and there are 109 failing pairs with m ≤ 20000.
