# Ge–Liu equiangular-line bound: exact checks

Candidate: Ge and Liu, *New bounds for equiangular lines and Balla’s conjecture*, arXiv:2606.29392v1, Conjecture 5.1 and Section 5, Question 2.

Run `python3 check.py` with Python 3 and SymPy. The recorded output is `check-output.txt`; success ends with `ALL_CHECKS_OK`. Every check uses integers, rational numbers or exact symbolic algebra.

The source fidelity checks reproduce the exceptional bounds 169 and 135 and its cycle excess identity. The geometry check constructs all points and isotropic lines of the symplectic generalized quadrangle over F₂, verifies the incidence identity and adjacency characteristic polynomial, and obtains a positive semidefinite Gram matrix of rank 21 for 30 lines at angle arccos(1/5).

The family checks verify symbolic multiplicity and dimension identities, the exact linear excess and uniform polynomial inequalities for k ≥ 8, plus integer evaluations for t = 3,…,12 with k = 2ᵗ and q = k²/2. In particular, k = 8 gives N = 67650 and d = 50226, exceeding the proposed branches 41328 and 56503.

The general incidence spectrum is classical: Ábrego, Fernández-Merchant, Neubauer and Watkins, *Trace-minimal graphs and D-optimal weighing designs*, Lemma 4.7, <https://www.csun.edu/~ba70714/publications/tracemin.pdf>. The Seidel/Gram transformation is classical; see Neumaier, *Graph representations, two-distance sets, and equiangular lines*, LAA 114–115 (1989), Proposition 2.1, DOI 10.1016/0024-3795(89)90456-4; it is also used by the candidate paper.

The full q ≥ 32 incidence matrices are not enumerated here. Their existence and spectra are supported by the classical finite-geometry construction and the mathematical derivation in the preregistration draft. These experiments are not a Lean proof, and do not establish exhaustive literature novelty.
