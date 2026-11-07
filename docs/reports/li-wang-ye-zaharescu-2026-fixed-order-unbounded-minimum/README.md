# Fixed-order minimum entries and entrywise cubes

Source: Xinyi Li, Shuochen Wang, Zimeng Ye and Alexandru Zaharescu,
*Simultaneous Unimodularity of Positive Integer Matrices and Their Entrywise
Cubes*, [arXiv:2610.14546v1](https://arxiv.org/abs/2610.14546v1), Section 9.
The question asks whether, for each fixed order at least ten, both determinants
can equal one while every entry exceeds an arbitrary prescribed bound.

Run `python3 check.py` with Python 3 and SymPy. `output.txt` is its output.
All arithmetic is integral, rational, or symbolic; no floating point is used.
The checker is independently written and does not import the paper's checker.

The source fidelity checks recompute the twenty signed seed cofactors with
fraction-free elimination and independent rational Gaussian elimination. They
verify the primitive-vector, Bezout, pairwise-gcd, homogeneous-point and
exceptional-prime certificates. The reduced real Jacobian minor is
1123741091386468277532. At the prime three the check uses the residual and
derivative valuations; it does not require a unit Jacobian minor modulo three.

The proposed argument's signed cofactor identity passes 3524 exact cases:
all 1-by-2 and 2-by-3 frames with entries in {1,2}, all relation coefficients
in {-1,1,2}, every insertion position, both ordinary and cubed frames, and
twenty additional tests on the source's 9-by-10 frame. The order-ten cofactor
homogeneity laws pass 60 cases at scales 2, 10 and 100. A symbolic identity
checks the scalar step used to exclude singularity after row replacement.
These tests include degenerate frames; they test algebraic identities, not
existence of a positive generic solution of the determinant equations.

**Missing required evidence:** no positive order-ten matrix with both
integer determinants equal to one has been computed for K = 2, 10 or 100.
The program prints each missing witness explicitly. Its final `ALL_CHECKS_OK`
certifies only the finite checks listed above. The source itself gives no
concrete order-ten solution; it proves existence by a circle-method theorem.
Neither the infinite existence claim nor that analytic theorem is verified by
this program. The preregistration draft requires revision before filing under
the requested numerical-witness requirement.

To validate supplied witnesses, run `python3 check.py witnesses.json`, with a
JSON object mapping the keys "2", "10" and "100" to integer 10-by-10 matrices.
The checker then verifies every strict entry bound and both determinants by
two independent exact algorithms, and fails on any invalid certificate.

The numerical source data come from Sections 2-4 of the cited TeX source.
The analytic counting theorem and the geometric ingredients require legitimate
proved Lean dependencies; no Lean proof or axiomatic substitute is supplied.
