# Equivariant resultant sign: exact finite checks

Source: Sonagnon Julien Owolabi, Ollé Grégroire Kam, Ibrahim Nonkane and Joel Tossa,
*Resultant of an equivariant polynomial system with respect to a direct product of r symmetric groups*,
arXiv:2609.26127v1, Remarks 3.8 and 5.5 and Theorem 5.4.

Run `python3 check.py`. Python 3 standard library only; all coefficients and determinants are
integers. `output.txt` contains the output. A failed assertion exits nonzero; success ends with
`ALL_CHECKS_OK`.

The fidelity check reproduces five exact coefficient specializations of the binary factor
in the source's Section 4 example with degree 3, n=4 and p=q=2. Two specializations vanish;
the other three give -3, -1863 and -31899.

The sign check specializes F_i=X_i^d for every ordered positive block composition of
2 <= n <= 5 and 1 <= d <= 3, and also (2,2,2) at d=1,2. This includes all six parameter
choices in Remark 3.8. It constructs divided differences by exact polynomial division
before identifying variables in each canonical partition piece. It independently compares
the resulting polynomials with the complete homogeneous sums, computes every ordered
partition resultant using the Macaulay matrices of Section 2.2, and checks the full unsigned
product and the blockwise coefficient-degree identity for the exponents. Leading constants
are checked separately; no degree-zero resultant is formed.

Sparse triangular matrices are checked entry by entry and their determinants are products
of their diagonal entries. Nontriangular determinants use fraction-free Bareiss elimination;
triangular matrices of dimension at most 60 are also cross-checked by Bareiss. Empty
extraneous matrices have determinant 1. Every Macaulay denominator is checked nonzero,
and every integer division is checked exact. A deliberate swap of the two linear
polynomials gives -1, checking sensitivity to the sign convention.

These are bounded exact computations, not a proof of the universal conjecture or a Lean
settlement. The classical coordinate-power specialization is due to Busé–Karasoulou,
arXiv:1407.2799, Section 3.3; no new technique is claimed here.
