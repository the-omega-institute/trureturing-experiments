# OEIS A367191: the conjectured formulas for the largest k with k² + kn + n² a square

A367191(n) is the largest positive k such that k² + kn + n² is a perfect square (0 if none). J.-F. Alcover conjectures formulas by n mod 4. They follow from the factorization u·v = 3n², where u = 2m − (2k+n) and v = 2m + (2k+n). Then k = (v − u − 2n)/4, which decreases as u grows, so the largest k comes from the least admissible u:
- u = 1 for n odd;
- u = 2 for n ≡ 2 (mod 4);
- u = 4 for n ≡ 0 (mod 4).

`check.py` (Python 3) does two things:
- for n ≤ N1 it compares a brute-force a(n) with the published data and the formulas;
- for n ≤ N2 it checks the formula's witness and that no smaller u is admissible.

Run `python3 check.py 150 100000`. Exit 0 and `RESULT PASS` mean every check holds.
