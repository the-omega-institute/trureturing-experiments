# Coherent-leakage minimizer checks

Source: Hassan Nasreddine, *Sharp Entropy Cost of Coherent Leakage Across a Block Decomposition*, arXiv:2605.11024v2, Remark IV.10 and Corollary IV.11. The source constructs a one-channel minimizer and leaves necessity open.

Run `python3 check.py` with SymPy installed. The stored output was produced with Python 3.14.8 and SymPy 1.14.0. All evidence uses integers, rational numbers and exact symbolic algebra; no floating-point comparisons occur.

The fidelity check reproduces the source's pure two-level spectrum and exact binary-entropy identity. The symbolic checks verify the determinant and antiderivative derivative identities, strict protected-mass monotonicity and the merging numerator decomposition used by the proposed proof. The rational grid checks the merging inequality and its exact equality condition for every enumerated tuple with one to three channels. Matrix certificates verify equation (14), its constraints and its spectra for dimensions one to three, including the singular endpoint and the saturated protected floor.

The finite grids check the proposed rigidity mechanism. They do not prove the integral identity, data processing, or the full minimizer classification over all complex matrices. The planned positive settlement uses a general analytic argument; these checks are supporting evidence.

The scalar enumeration uses protected masses 1/2 and 1, leakage masses 0, 1/4 and 1/2, coherence fractions 0, 1/2 and 1, residual leakage 0 and 1/4, and resolvent parameters 1/3, 1 and 2. The matrix enumeration uses both dimensions 1 through 3, floor 1/(2d_P+2), leakage 1/4, 3/8 or 1-d_P a_0, coherence fractions 1/4, 1/2 or 1, rational unit vectors, and phases 1, i or (3+4i)/5.
