# Isoclinism and faithful complex representation dimension

This is an exact sanity check for Section 5, Question (2) of Gurleen Kaur,
Amit Kulshrestha and Ayush Udeep, *On the representation dimension of finite
p-groups*, [arXiv:2608.18814v1](https://arxiv.org/abs/2608.18814v1).
The question concerns all finite groups of equal order with isomorphic centers
and an isoclinism, not only the groups of order p^6 in the motivating table.

Run `python3 check.py` with Python 3.12 or newer and SymPy 1.14.0.
The program uses integer polynomial arithmetic modulo cyclotomic polynomials.
It constructs actual matrices for every irreducible representation in the
standard complete list of the 28 group isomorphism types of orders 1 through 15,
eight additional direct products of orders 16, 24 and 32, and a second
presentation of S3 x C2 used for source fidelity. The completeness of the
28-type group list is the standard small-order classification input; the
program does not enumerate every possible multiplication table to prove that
classification. It verifies each multiplication table's group laws and each
matrix homomorphism, character orthogonality and the degree-square sum.
Thus its irreducible lists are complete for each supplied group. Kernel-mask
dynamic programming computes the exact least dimension of a faithful direct
sum. The source's positive-dimension convention assigns dimension 1 to the
trivial group.

For every same-order pair, the program searches for a center isomorphism and
for a central-quotient isomorphism inducing a derived-subgroup isomorphism on
commutators. It tests the resulting identity on all 45 eligible unordered pairs,
including self-pairs, among the 37 supplied presentations. Nonisomorphic pairs
include D8 and Q8, D12 and Dic12, and their listed central direct factors.
The non-nilpotent D12/Dic12 case is included.

Source fidelity reproduces the Section 3 example
delta(S3 x C2) = 2 < delta(S3) + delta(C2) = 3, as well as the Section 1
statement that C2 x C2 has no faithful irreducible representation.
The output is in `check.out`; its final line is `ALL_CHECKS_OK`.

These finite checks support the preregistered general proof route. They prove
neither the universal assertion nor the completeness of the prior-art screen.
No Lean proof of the assertion or its negation is included.
