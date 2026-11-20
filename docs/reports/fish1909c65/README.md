# Fish et al., Conjecture 6.5

Run `python3 check.py` with Python 3.12 or newer and SymPy 1.14.0. The program uses fractions and exact arithmetic in real quadratic fields. Every comparison with a square root reduces to a comparison of rationals.

The fidelity check reproduces the Chebyshev-norm configurations in Example 4.5(1) for sizes 4 through 8. The counterexample check covers sizes 4 through 12 in the norm given by the maximum of the Euclidean norm and twice the absolute vertical coordinate. It verifies distance classes, the slope and height bounds for the circle argument, noncollinearity and the absence of arithmetic-progression projections. The projection check tests affine interpolation on the first four points against all possible distinct ranks in a progression of the specified size. The preregistration contains the argument for unbounded size. [arXiv:1909.08769v2](https://arxiv.org/abs/1909.08769v2)
