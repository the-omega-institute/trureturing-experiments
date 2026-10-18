# OEIS A246056: residues mod 3 through a Legendre generating-function identity

`check.py` computes A246056 from its generating function, compares with the entry's b-file (`b246056.txt`, copied from oeis.org), checks the exact identity Sum_n x^n (1+x)^(-2n-1) Sum_k C(n,k)^2 (-x)^k = 1/sqrt(1+4x^2), the termwise congruence of the generating function to that series mod 3, and the residues a(n) = 1 (mod 3) iff n = 2m with m free of the base-3 digit 2 (else 0): the settlement of the open-problem preregistration in `the-omega-institute/trureturing`.

Run: `python3 check.py 300 b246056.txt` (exit 0, `RESULT PASS`; output in `out.txt`).
