# Open-ended QAOA positivity on P4

Exact arithmetic check of Conjecture 10 in Inbar Chefer, Uri Shaham and Adi Makmal,
“Measurements Number Scaling in the Quantum Approximate Optimization Algorithm for MaxCut:
A Statistical Analysis”, [arXiv:2607.03340v1](https://arxiv.org/abs/2607.03340v1).

Run `python3 check.py` with Python 3 and SymPy 1.14.0. The stored output is in
[check.out](check.out); a successful run ends with `ALL_CHECKS_OK`.
All evidence uses exact algebraic arithmetic, with no floating point.

The graph has vertices 0,1,2,3 and edges 01,02,23. The root is 01 and the
open-ended two-edge continuation is 0,2,3. Use the source convention
C = −Σ C_e, phase exp(+i γ Σ C_e), mixer exp(−i β Σ X_i), and basis |x0 x1 x2 x3>.
At (γ1,γ2,β1,β2) = (π/2,−π,−3π/8,3π/8), the computed normalized state is
(|0000> − i|0011> − i|1100> + |1111>)/2. The root cut expectation is exactly zero.

All four angles are also perturbed by
(2 atan(1/100),2 atan(1/97),atan(1/101),atan(1/103)). The root expectation becomes
4532808865055386329699672853/6485560649008405461294892502500 > 0,
and all sixteen amplitudes are nonzero. The conjecture quantifies over every
parameter vector; it imposes no generic-angle restriction. The exceptional zero
is a small counterexample to that literal quantifier. These tests do not refute
an almost-everywhere formulation or the paper's repeatable-geodesic lemma.

Fidelity checks reproduce source Lemma 7 on these two depth-two states: global
bit-flip invariance and all single-qubit Z expectations zero. A separate exact
depth-one P4 state checks distant factorization: vertices 1 and 3 have distance
3 > 2p, and their ZZ expectation vanishes. A zero-angle check gives every
local cut expectation 1/2. These are finite instance checks, not proofs of Lemma 7.

The program constructs each full 16-by-16 mixer tensor matrix and each diagonal
phase layer. SymPy's exact algebraic number field canonicalizes the radicals.
Rational-tangent perturbations keep all computations exact; at four qubits the
four single-qubit normalization denominators multiply to an integer square.
The result is experimental evidence for a later Lean refutation, not a Lean proof.
