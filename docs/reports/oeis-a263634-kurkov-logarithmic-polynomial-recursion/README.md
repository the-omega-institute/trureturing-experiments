# Kurkov's logarithmic-polynomial recursion (OEIS A263634)

Source: [OEIS A263634](https://oeis.org/A263634), revision 166, modified
2025-11-05; the FORMULA conjecture is credited to Mikhail Kurkov, 2025-03-30.
The coefficient ring is Q[x1,x2,...]. The conjecture identifies R(n,1), where
R(1,k)=x_k and R(n,k)=R(n-1,k+1) minus the convolution with weights
binomial(n-2,j-1), with n! times the coefficient of t^n in
log(1 + sum of x_i t^i/i!).

Run with Python 3 and SymPy:

    python3 docs/reports/oeis-a263634-kurkov-logarithmic-polynomial-recursion/check.py

The saved output is `check.out`. All arithmetic is exact over Q. The program:

- Reproduces the six polynomials and coefficient rows printed in EXAMPLE.
- Checks the full symbolic identity for every 1 <= n <= 8, with all variables
  symbolic; the term counts are 1, 2, 3, 5, 7, 11, 15, 22.
- Checks the source's row-sum and partition-number assertions through order 8.
- Independently checks R(n,k) against (n-1)! times the coefficient of t^(n-1)
  in (D^k F)/F for all 36 positive pairs with n+k-1 <= 8.

These finite checks are evidence for source fidelity and the proposed proof
route. They do not prove the universal identity. The program prints
`ALL_CHECKS_OK` only after every assertion succeeds.
