# OEIS A356334: solutions of x^(n+1) + y^(n+1) = (x+y)^n

A356334(n) counts the solutions 0 ≤ x ≤ y of x^(n+1) + y^(n+1) = (x+y)^n. C. Wu conjectures a(n) = 3 except for n = 0 and n = 2.

**Proof.** Write x = da and y = db with gcd(a, b) = 1. Then d(a^(n+1) + b^(n+1)) = (a+b)^n.
- **n odd.** The gcd of a + b and a^(n+1) + b^(n+1) divides 2, so a^(n+1) + b^(n+1) is a power of 2.
- **n even.** R = (a^(n+1) + b^(n+1))/(a+b) is at most n + 1 (lifting the exponent), but at least b^n/2.

Only (0, 0), (0, 1) and (2^(n−1), 2^(n−1)) remain, plus (1, 2) at n = 2.

`check.py` (Python 3) finds all solutions for n ≤ NMAX by brute force over the range the proof allows. It also checks the inequality R > n + 1 for even 4 ≤ n ≤ 60 and coprime 1 ≤ a < b ≤ 60.

Run `python3 check.py 14`. Exit 0 and `RESULT PASS` mean every check holds (about 55 s).
