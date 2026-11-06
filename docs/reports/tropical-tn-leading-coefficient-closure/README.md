# Tropical TN leading-coefficient closure

The exact check supports the proposed answer to the first further-study
question in Section 6 of Saeidi Dahaki and Manjegani,
*Spectral Equality for Tropical Totally Nonnegative Matrices: Dominant
Cancellation and Newton Polygons*, arXiv:2610.06887v1:
the admissible leading-coefficient space contains a full open box and hence
has zero vanishing ideal in its finite-entry coordinates.

Run `python3 check.py`. The program uses only Python's standard library,
integers, and exact rational arithmetic. Its output is in `check.out`.

The fidelity checks reproduce Example 5.3's dominant coefficient 62 and
Example 5.17's coefficient polynomials, determinant expansion, and exceptional
lifts at three rational epsilon values. The tied 2 by 2 case includes both
cancellation and strictly positive dominant determinants at the same
valuation, so a tropical tie does not force equality of leading products.

For dimensions 1 through 3 the program exhausts all matrices with finite
diagonal in {-1,0,1} and other entries in {absent,-1,0,1}. It also exhausts
all 4 by 4 supports with zero finite weights and finite diagonal. It checks
ordered maximizing matchings, strict integer rearrangement gaps, and the
exact coefficient bound on the entire proposed box. It separately evaluates
all 512 corners of the 3 by 3 zero-weight box.

The box uses b = 1 + 3^n n! and center c_ij = b^(ij), with indices starting
at 1 and each coordinate between half and three halves of its center.
These are finite exact tests of a general proof route. They do not establish
the universal closure theorem, and no Lean proof is supplied here.
