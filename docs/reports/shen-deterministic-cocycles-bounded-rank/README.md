# Deterministic cocycles with bounded contraction radicals

`check.py` uses Python 3 standard-library integer arithmetic and `Fraction`.
Run `python3 check.py` to reproduce `output.txt`; successful completion ends
with `ALL_CHECKS_OK`.

The fidelity check exhaustively reproduces all six rank-distribution rows of
the finite-rank census in Ce Shen, *Maximal Total Quantum Dimension at Bounded
Rank in Dijkgraaf-Witten Theories*, arXiv:2609.26914v1, Appendix A.2, and checks
their averages against Theorem 3.8 using exact fractions.

The proposed construction on K = F_(p^n) is the prime-field trace of the 3 by 3
Moore determinant with rows of exponents 1, p, p². The program implements
finite fields as polynomial quotients, certifies the selected moduli with
Rabin's exact irreducibility criterion, and checks trace nondegeneracy,
Frobenius linearity, alternation, contraction matrices, radical sizes, the
trace-adjoint expression, and the signed degree-p⁴ contraction polynomial.
It checks every ordered pair in each of 19 fields: p=2 with 1≤n≤8, p=3 with
1≤n≤5, and p=5,7 with 1≤n≤3, totaling 290,183 pairs. Contractions are evaluated
at every vector of the third-slot polynomial basis and compared with the
matrix evaluation. Coordinate-cochain alternation is checked on every basis
triple, and on every field triple when the field has at most 16 elements.
The five-term additive cocycle identity is checked on every quadruple when
the field has at most nine elements. Sparse bivariate polynomial arithmetic
also checks the signed expansion at p=2,3,5,7,11,13.

All tested nonzero-flux radicals have at most p⁴ elements. All computed rank
sums satisfy (p+1)p^n−p ≤ r < (p⁴+1)p^n. In the tested range n≥5 the forms
have zero full radical. The p=3,n=3 form is zero, as its computed rank 729
records; the construction claims nondegeneracy only for n≥5.

These are finite checks, not a proof for all dimensions or all primes. The
modulus search in this check is exhaustive and is not an efficient field
construction. The proposed efficiency argument instead uses Victor Shoup,
*New Algorithms for Finding Irreducible Polynomials over Finite Fields*,
Mathematics of Computation 54 (1990), 435–447,
DOI: 10.1090/S0025-5718-1990-0993933-0, for fixed p.

The binary degree-16 radical bound has prior art: Mesnager–Kim–Jo,
arXiv:1811.11280, Section 7; Gao–Kan–Li–Xu–Wang, arXiv:2309.11229,
Theorem 6. It is not claimed as new. The candidate concerns the signed
all-prime family and its application to Shen's Section 5 first outlook
question. Sources supply mathematical statements; this check program is
independently written.
