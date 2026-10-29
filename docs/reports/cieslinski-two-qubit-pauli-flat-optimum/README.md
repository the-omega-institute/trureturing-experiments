# Two-qubit Pauli-flat optimum

Source: Paweł Cieśliński, *Pauli Flat Quantum States Mimicking Maximal Magic*, [arXiv:2610.07172v1](https://arxiv.org/abs/2610.07172v1), Definition 1, equations (6)–(7), Section V optimization (15), and Table 1.

Run `python3 check.py`; Python 3's standard library suffices. The expected output is in [output.txt](output.txt), ending in `ALL_CHECKS_OK`.

The fidelity check reconstructs the source's displayed two-qubit state at gamma = 1/3. It checks Hermiticity, trace one, all 15 principal minors, all 15 Pauli expectations, its characteristic polynomial, and purity 2/3. Table 1 prints 0.3334 and 0.6667; these printed approximations are not used as exact feasible parameters.

The independent check verifies the six Mermin-square operator identities and all 24 associated positive rank-one projectors using Gaussian integers. It enumerates every one of the 32,768 arbitrary sign assignments, finds a parity-mismatched context, and checks Tr(T) = 4 and Tr(T H) = -12. Thus tau = T/4 has expectation (1 - 3 gamma)/4 in rho = (I + gamma H)/4; positivity forces gamma <= 1/3. It also checks Tr(H^2) = 60 for every assignment and an explicit Gram matrix attaining gamma = 1/3 and purity 2/3. The real parameter gamma is handled by these exact affine coefficients, not by a sampled grid.

All arithmetic uses integers, Gaussian integer pairs, and fractions. No floating point or numerical semidefinite solver is used. The result is an exact computational check and a mathematical proof outline; it is not a Lean proof. It addresses N = 2 only.
