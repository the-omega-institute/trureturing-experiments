# Clarke Conjecture 8.7: exact finite checks

Source: Oliver Clarke, *The equivariant Ehrhart theory of polytopes with
order-two symmetries: Part II*, [arXiv:2610.17988v1](https://arxiv.org/abs/2610.17988v1),
Conjecture 8.7 and Sections 2.1–2.2. The conjecture asks whether effectiveness
of the equivariant numerator persists under every positive integral dilation.

Run `python3 check.py > check.out` with Python 3 and SymPy. The stored output
was generated with SymPy 1.14.0. All evidence uses integers and exact rational
linear algebra; no floating-point approximation contributes to a check.

The representation ring of the group of order two is represented by pairs
of trivial/sign multiplicities, with tensor multiplication
`(a,b)*(c,d) = (ac+bd,ad+bc)`. The checker verifies 125 exact polynomial
identities `D(t) K_m(t) = D(t^m)` for `1 <= p,q <= 5` and odd `1 <= m <= 9`,
including coefficientwise effectiveness of the geometric-sum kernel.

Geometric fidelity checks reproduce the source's orthogonally and diagonally
reflected squares (Examples 2.3–2.4), pyramid, bipyramid, and effective polytope
with an obstructing face (Examples 5.1–5.3). The last example is four-dimensional;
its integer fibers are counted directly from the tetrahedron inequalities.
Twenty-seven additional random symmetric lattice polytopes have dimensions
one through three, with seed 261087. Their vertices and actions are printed
in the output so every case can be reproduced.

Facets are obtained from exact rational nullspaces, and lattice-point counts
use exact integer intervals along coordinate fibers. An independent exhaustive
point enumeration cross-checks every polytope at dilation parameters 0–3
(128 comparisons). Counts through parameter 56 produce numerator coefficients
for dilations 1–7 through degree 8 and verify 128 odd-dilation extraction
identities. Among the 32 geometric examples, 24 have effective base coefficients
through degree 56 and eight are negative controls. All 224 dilation cases pass:
checked-effective bases retain effectiveness in the tested window, and every
even dilation is effective in that window. Complete coefficients appear in
`check.out`, whose final line is `ALL_CHECKS_OK`.

These are finite checks, not a proof of the conjecture, of polynomiality of an
arbitrary numerator, or of effectiveness beyond the stated coefficient windows.
The proposed settlement uses a general semiring argument for odd dilations
and the source's no-obstructing-face theorem for even dilations. No Lean proof
is part of this report.
