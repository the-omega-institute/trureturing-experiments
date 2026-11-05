# Alstad–Duff–Katzman image-equation checks

This report checks the set-theoretic reading of Example 4.6 after Equation (38)
in Colin Alstad, Timothy Duff and Mordechai Katzman, *Segre-Determinantal Loci
and the Image Variety for Three Flatland Cameras*, arXiv:2606.31965v1.
The target concerns two variable cameras from projective 3-space to projective
2-space over the complex numbers, for every number of correspondences greater
than eight. The source explicitly includes set-theoretic equations among its
requested descriptions. Ideal generation is outside this report's claim.

Run `python3 check.py` with Python 3 and SymPy (tested with SymPy 1.14.0).
All arithmetic is integer, rational or symbolic exact algebra. The deterministic
seed is 260631965. Stored stdout is in `output.txt`; a successful run ends with
`ALL_CHECKS_OK` and exits 0.

The program reproduces the source's fundamental block determinant, its
vectorized epipolar constraint and eight-point cofactor construction. It tests
all eight-row cofactor determinants and nine-row Segre minors on actual camera
tuples with 8, 9 and 10 correspondences. A negative control with kernel spanned
by the identity matrix is rejected. Boundary controls include a rank-eight
Segre matrix whose kernel consists only of multiples of a rank-one matrix,
and a rank-one Segre matrix with identically zero eight-row cofactors.

For the proposed additional argument it verifies the three polynomial
rank-one-to-rank-two deformation formulas symbolically and on all 88 compatible
pairs in a 13-point rational projective grid. It checks reconstruction formulas
for all 61 pairs in the grid satisfying the normalized rank-two epipolar
constraint, including both epipoles. Each curve can use the same parameter
across any finite tuple. The rank-one deformation formulas are also checked
with independent symbolic parameters.

These checks do not prove the statement for all numbers of correspondences,
do not certify Zariski-closure membership, and are not a Lean proof. The
preregistration's proposed settlement uses a general argument; the bounded
checks are diagnostics, not its main mathematical content.
