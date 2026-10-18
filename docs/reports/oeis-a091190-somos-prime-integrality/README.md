# OEIS A091190 / A091200: x·A(x)^p = B(x·A(x^p)), B(x) = x/(1 − px), is integral exactly for prime p

Comment in both entries: "More generally, given A(x) satisfies x*A(x)^p = B(x*A(x^p)) where B(x) = x/(1-p*x), then it appears that A(x) is an integer series only when p is prime."

## Route

The equation is A(x)^p = A(x^p)/(1 − p·x·A(x^p)) with A(0) = 1, and the coefficient of xⁿ gives p·a_n = [xⁿ](A(x^p)/(1 − p·x·A(x^p)) − A_{<n}(x)^p).
- **p prime.** The bracket is ≡ 0 (mod p), because A_{<n}(x)^p ≡ A_{<n}(x^p) (mod p). Hence A is integral.
- **p composite, least prime factor q.**
  - For n < p, A ≡ (1 − px)^{−1/p}, so a_n = Π_{i<n}(1 + ip)/n!.
  - Every factor 1 + ip is ≡ 1 (mod q), while q! contains q. So a_q is not an integer.

## Program

`check.py` (standard library only; reads the two entry texts saved here). Run `python3 check.py`; it prints `RESULT PASS` when every check below holds.
- **Primes p ≤ 29.** A is integral to order 14 (p < 10) or 8.
- **Composite p ≤ 30.** The first non-integer coefficient is a_q, and it equals the closed form.
- **The entries.** The first 14 terms of A091190 (p = 3) and A091200 (p = 5) match.

Reading: `run.log` (exit 0, `RESULT PASS`).
