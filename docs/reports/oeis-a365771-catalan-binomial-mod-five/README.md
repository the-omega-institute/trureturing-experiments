# OEIS A365771: Hanna's three congruences modulo 5

a(n) = C(2n+1, n)/(2n+1)·C(3n−1, n). With A033042 the numbers whose base-5 digits are 0 or 1, P. D. Hanna conjectures:
- a(5m) ≡ 4 (mod 5) for m ∈ A033042, m > 0;
- a(5m+1) ≡ 2 (mod 5) for m ∈ A033042, m > 0;
- a(n) ≡ 0 (mod 5) otherwise (n > 0).

For n ≥ 1, a(n) = 2M(n)/(3(n+1)) with M(n) = (3n)!/(n!)³. Then:
- **Kummer.** v₅(M(n)) is the total carry of n + n + n in base 5, which exceeds v₅(n+1) as soon as n has a digit ≥ 2.
- **Lucas.** M(n) ≡ 1 (mod 5) when all digits are 0 or 1, so a(n) ≡ 2/(3(n+1)) ∈ {4, 2}.

`check.py` (Python 3) checks:
- the identity and the published data;
- the valuation comparison and M(n) mod 5;
- the three conjectures for n ≤ N.

Run `python3 check.py 6000`. Exit 0 and `RESULT PASS` mean every check holds (about 20 s).
