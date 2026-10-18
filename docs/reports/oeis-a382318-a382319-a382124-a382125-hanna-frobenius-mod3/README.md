# OEIS A382318/A382319 and A382124/A382125: residues mod 3 through the Frobenius congruence

`check.py` computes the four sequences from their definitions, compares them with the entries' b-files (copied from oeis.org), and checks the steps of the open-problem preregistrations in `the-omega-institute/trureturing`:

- **A382318/A382319.** sum_{n>=1} A(x^n) = B, A = x + B^3 and B^3 = B(x^3) (mod 3), hence a(3n) = b(n) and a(3n+1) = a(3n+2) = 0 (mod 3).
- **A382124/A382125.** sigma(n)sigma(2n)/3 is an integer satisfying the Dold congruences (checked for p < 40), C = A382124 is integral, D = C^3 and C^3 = C(x^3) (mod 3), hence d(3n) = c(n) and d(3n+1) = d(3n+2) = 0 (mod 3).

Run: `python3 check.py 300` (exit 0, `RESULT PASS`; output in `out.txt`).
