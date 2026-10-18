# OEIS A356508: residues mod 4 of the solution of 2 = Product (1 + x^n A)(1 + x^(n-1)/A)

`check.py` solves 1/A = 2/P - 1 (P the product without its factor 1 + 1/A) by its division-free recursion, checks the defining product, compares with the entry's b-file (`b356508.txt`, copied from oeis.org) and checks the steps of the open-problem preregistration in `the-omega-institute/trureturing`: A = 1 (mod 2), each factor (1 + x^k A)(1 + x^k/A) = (1 + x^k)^2 (mod 4), 1/A = 1 + 2 sum p(n) x^(2n) (mod 4), and the residues of a(n).

Run: `python3 check.py 120 b356508.txt` (exit 0, `RESULT PASS`; output in `out.txt`).
