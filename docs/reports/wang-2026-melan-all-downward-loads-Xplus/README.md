# Wang's Melan equation: exact preparatory checks

Target: Jinxiang Wang, *On the well-posedness and efficient approximation for
the classical Melan equation in suspension bridges*, arXiv:2609.08912v1,
Section 5's heavy-load question, interpreted in the source's nonnegative-integral
deflection class X_+. Source: https://arxiv.org/abs/2609.08912v1.

Run `python3 check.py` with Python 3.12 or newer and SymPy 1.14.0. The stored
output is `check.out`; successful execution ends with `ALL_CHECKS_OK`.
The program uses integers, rational numbers, exact exponential expressions,
symbolic differentiation and integration, and rational-polynomial Sturm counts.
It uses no floating-point comparisons.

The fidelity check reproduces the source's Green formulas (2.2), (2.3), (2.5),
including their derivative jump and all four hinged boundary traces. It
reproduces the double integral appearing in (2.10) and the constant-load case
of the linear solution formula (2.6), for three exact parameter choices.

The proposed continuum mechanism is checked on Dirichlet finite-difference
matrices. Both resolvent identities and their coordinate signs are checked in
144 cases with one through six interior nodes. Nine heavy-load models use
two, three, or five nodes and uniform, left-span, or central loads of amplitudes
100, 10000, or 1000000. For each, the scalar balance has exactly one positive
root by a Sturm count, and exact rational bisection encloses it after 32 steps.
The reconstructed frozen linear equation has zero exact residual. A rational
Taylor enclosure for exp(1) proves the source's (H) threshold is below 3 for
a=b=c=L=1; all tested discrete load masses exceed 3.

These checks support the proposed signs, identities, and uniqueness mechanism.
They do not prove the continuum theorem, Sobolev reconstruction, convergence
of the discretizations, or convergence of the source's Picard iteration.
The positive settlement must use a general continuum argument; these finite
checks are preparatory evidence only.

The program is original code derived from the displayed mathematical formulas.
SymPy is distributed under BSD-3-Clause; its implementation is not copied.
