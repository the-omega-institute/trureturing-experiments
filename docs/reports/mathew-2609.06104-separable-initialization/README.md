# Separable initialization for Mathew's complete-graph quantum battery

Source: Sheryl Mathew, *Battery Locality Is Necessary in Noncommuting Quantum Charging Bounds*, [arXiv:2609.06104v1](https://arxiv.org/abs/2609.06104v1), equations (11)–(23) and (106)–(113), and the concluding initialization question.

`check.py` uses exact SymPy coefficients and a sparse Pauli algebra for the literal edge-qubit battery and explicit Jordan–Wigner auxiliary Majoranas. It uses no floating-point arithmetic. A Pauli word is stored as a coefficient times X^x Z^z, with multiplication sign determined by the binary dot product z·x′. The trace is the identity coefficient times the Hilbert-space dimension.

For m = 1, 2, 3, 4 it reproduces the source identities H_B² = 2m I and C² = m I. It checks the Clifford relations, i[C,H_B⊗I] = 2√2 D, {C,D} = 0, [C,D²] = 0, and the full second-moment expansion. It constructs nonzero commuting spectral projectors for battery states in the ground sector and opposite auxiliary paired eigenvalues, normalizes them to density matrices, and checks the proposed ground-sector product peak at t = π/(4√m). The projectors are Hermitian idempotents, so their normalized positive traces supply positive density matrices. The resulting covariance inner products are 0, −2, −2, −4; the peaks are 2√2, 4√2, 10√2/3, 4√2. These finite checks support a proposed general proof and do not prove the universal separable inequality.

Run `python3 check.py`; its stored output is `check.out` and its last line is ALL_CHECKS_OK.

`numeric_check.py` is separate exploratory floating-point data, using NumPy. For each m = 1, 2 it samples 1,024 independent random pure battery–auxiliary products and maximizes the sinusoidal power over all times analytically. It also performs 20 alternating largest-eigenvector updates from 32 random product starts, updating the time phase each sweep. Seeds are 260906105 and 260906106. The sampled time maxima are 2.825662280578 and 1.272981272523. The alternating searches reach 2.828427124746 and 5.656854249492, matching the proposed bounds within floating-point rounding. The benchmark values are 4 and 5.656854249492. These searches do not certify global maxima or prove a bound. Their stored output is `numeric_check.out`.

Run the exploratory search with `OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 numeric_check.py`.

Runtime used for the stored outputs: Python 3.14.8, SymPy 1.14.0, NumPy 2.4.2.
