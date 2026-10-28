# KD-real mixed-state SLD counterexample

Run `python3 check.py` with SymPy installed. The recorded run uses SymPy 1.14.0.
The script uses integer modular phases and exact symbolic algebra; it never uses
floating point to decide an identity or inequality.

The source is Spriet–Langrenez–Brummelhuis–De Bièvre,
[arXiv:2511.01996v2](https://arxiv.org/abs/2511.01996v2), Section 6.3 and Appendix E.
The fidelity check reproduces Proposition 6.5 on a pure qubit: QFI is 4, the SLD
has a non-real KD symbol, and its conditional expectation is nonzero and
anti-Hermitian. It also checks Proposition 6.3 on the mixed family: the two basis
measurements have zero classical Fisher information.

For `G = F_3^2`, the script checks the full symbolic interval `0 < t < 1` for
`rho_t = (I+t R)/9`. It checks a rank-three projection decomposition, exact lower
bounds on both density eigenvalues, strictly positive KD entries, noncommutation,
the SLD equation, reality of the SLD symbol, and `QFI = 3*t^2/(4+t)`.
At `t=1/2`, the density eigenvalues are `1/6` and `1/12`, the KD entries are `1/54`
and `1/108`, and QFI is `1/6`. The inequalities on the interval follow from the
printed rational formulas and the explicit positive factors `t`, `1-t`, and `4+t`;
they are not inferred from numerical samples.

For each of `p = 3, 5, 7`, all monomial Weyl products, adjoints, opposite-index KD
phase identities, and trace-orthogonality support conditions are checked with
integers modulo `p`. These finite checks do not prove the all-prime family. The
arbitrary-prime argument uses the antisymmetric Weyl algebra and invertibility
of its Jordan map for a faithful density matrix. Its Lean proof is outstanding.

`output.txt` is the exact program's output. Success ends with `ALL_CHECKS_OK`.
This directory contains experimental evidence, not a kernel-checked settlement.
