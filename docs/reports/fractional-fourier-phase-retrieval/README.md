# Exact checks for the fractional Fourier phase-retrieval conjecture

The target is Conjecture 13 of A. S. Bandeira, J. Cahill, D. G. Mixon and A. A.
Nelson, *Saving phase: Injectivity and stability for phase retrieval*,
arXiv:1302.4618v3, Section 2.2. The defining transform is Ç. Candan, M. A. Kutay
and H. M. Ozaktas, *The discrete fractional Fourier transform*, IEEE TSP 48(5),
1329–1337, DOI 10.1109/78.839980, equations (16), (22) and Table IV.

Run `python3 check.py`; Python 3 and SymPy are required. `check.out` is the stored
output. The program uses exact arithmetic in rational algebraic fields and
symbolic identities. It uses no floating-point evidence or numerical eigensolver.

The source-fidelity check reproduces the successful HMW test for the source's
3-by-8 matrix in equation (8): its Hermitian analysis map has rank 8 and its
one-dimensional nullspace has a normalized generator with determinant 9/8.

For dimension 4, the program supplies the entire real orthonormal Candan basis,
checks S eigenvectors, cyclic parity, within-parity eigenvalue ordering and the
ordinary DFT, and uses the correct mode labels 0, 1, 2, 4. The orthogonal unit
states (u_0 ± i u_4)/sqrt(2) have no zero coordinates and give exactly
(3/8, 1/8, 3/8, 1/8) in all four prescribed bases. Their pure-state difference
is a rank-two Hermitian nullspace element with eigenvalues +1 and -1.

The supplied search candidate's claim that dimension 3 remains injective is
false for the displayed source matrix. Its Hermitian analysis rank is 6 and
its nullspace contains a rank-two matrix. The real orthogonal unit states
(1/sqrt(2), 1/2, -1/2) and (1/sqrt(2), -1/2, 1/2) have equal intensities at
orders 0, 1/2, 1, 3/2: respectively (1/2,1/4,1/4), (1/3,1/3,1/3),
(1/6,5/12,5/12), (1/3,1/3,1/3). The program reports this source discrepancy
explicitly. The successful fidelity example is the separate matrix in equation
(8), not this incorrectly described fractional-transform example.

A symbolic real-coordinate identity checks the two-mode sign ambiguity. A
symbolic integer identity checks the affine-grid relative phase (-1)^p for
orders a_0+p/Q and mode gap 2Q. Finite index checks cover M=4 through 128 and
320 boundary/nearby cases for Q=1 through 32. Those finite checks do not prove
arbitrary-dimensional source-basis existence; that input is the Candan
construction. The general algebraic argument belongs in the preregistration
outline. These experiments are not a Lean proof or an exhaustive priority claim.

Primary sources:

- https://arxiv.org/pdf/1302.4618v3
- https://users.metu.edu.tr/ccandan/pub_dir/dis_frft_SP_TRAN00.pdf

`check.py` is an independently written check program. No source-paper file is
redistributed in this report directory.
