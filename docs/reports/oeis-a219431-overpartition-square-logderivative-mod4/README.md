# OEIS A219431: residues mod 4 of the logarithmic derivative of sum pbar(n^2) x^n

`check.py` computes overpartition numbers from their product generating function, checks pbar(m) = 2 d(m) (mod 4), computes G(x) = sum pbar(n^2) x^n (A219430) and the logarithmic derivative x G'/G (A219431), compares with the entry's b-file (`b219431.txt`, copied from oeis.org) and checks a(n) = 2n (mod 4), the settlement of the open-problem preregistration in `the-omega-institute/trureturing`.

Run: `python3 check.py 200 b219431.txt` (exit 0, `RESULT PASS`; output in `out.txt`).
