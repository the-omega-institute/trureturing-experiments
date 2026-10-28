# Semimonotone exact order and inverse Z-matrices

The source is B. P. Chauhan and D. Dubey, *On semimonotone matrices of exact
order two*, arXiv:2603.00639v1 (submitted 28 February 2026), Definition 3.1,
Example 3.1 and Remark 4.1, second bullet. The target is the general-order
question, not Conjecture 3.12 about exact order two.

Run `python3 docs/reports/semimonotone-inverse-z/check.py` from the repository
root with Python 3 and SymPy. The recorded output is `output.txt`. All tests use
integers, rational numbers and symbolic polynomial or rational identities.

For integers s,m >= 2 set L=3s^2, n=mL, d=3s(s-1), a=3s-1 and b=3s-2.
Partition the indices into m clusters of size L. The real symmetric matrix has
diagonal d, within-cluster off-diagonal entries -a and between-cluster entries
-b. The small-block maximal-coordinate certificate has margin
`d-(s-1)a=s-1>0`; every larger principal block has negative row sums since
`d-sb=-s<0`. These certificates give exact order k=n-s in both the ordinary
and strict senses.

With P the averaging projection on each cluster and Q the global averaging
projection, the matrix is `(L-1)(I-P)-(P-Q)-C Q`, where C=1+nb. Its proposed
inverse is `(I-P)/(L-1)-(P-Q)-Q/C`. The between-cluster entries of the inverse
are b/C>0, and the determinant formula is `(-1)^m (L-1)^(n-m) C`.

The program checks the parameter identities symbolically in the nested
projection algebra. Independently constructed full matrices at (s,m)=(2,2),
(2,3),(3,2) have exact inverse and determinant checks, projector identities,
and exhaustive principal-block occupancy types of every nonempty order.
The first case also has an independent exact characteristic-polynomial
check. The literal source matrices in Example 3.1 supply the fidelity check.
An additional 361 parameter pairs check the sign and cutoff identities.

The checks support a preregistration outline. They do not constitute a Lean
proof or establish literature novelty. No fixed-k unbounded-dimension claim
is made, and no conclusion about Conjecture 3.12 is asserted.
