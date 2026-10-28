# Vidal-tail equality classification: exact checks

The target is the unnumbered equality-state classification question following
Corollary 1 in Section III.1.3 of Xiong, Jin, Wang, Chen, Tao and Sze,
*Multiscale Schmidt-Spectrum Bounds for High-Dimensional Entanglement:
Geometric Measures and Schmidt-Number Witnesses*, arXiv:2609.21209v1.
Source: https://arxiv.org/html/2609.21209v1 (bound (14)).

Run `python3 check.py` with Python 3 and SymPy 1.14.0. All arithmetic uses
integers and exact algebraic expressions. The recorded standard output is
`output.txt`; success requires exit 0 and final line `ALL_CHECKS_OK`.

The check reproduces the source's isotropic-state fidelity identity for
symbolic fidelity in dimensions 2 through 5, including its dimension-four
benchmark. It verifies the proposed nonmaximal-reference mixed family for a
symbolic mixing parameter, its rank-two criterion and half-mixture spectrum,
and real and complex cutoff rotations over a finite grid of Schmidt spectra
and scaling parameters. Rectangular zero-padding is also checked.

For the explicit mixed family the exact fidelity and the displayed ensemble
cost coincide with the source's lower bound. Inferring the convex-roof value
uses that published bound. The program does not prove the lower bound,
convex-roof attainment, the universal classification, or absence of prior art.
Those remain obligations of the proposed mathematical settlement. No Lean
proof is included.
