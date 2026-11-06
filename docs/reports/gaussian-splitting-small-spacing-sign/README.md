# Gaussian splitting: the small-spacing coefficient

Run `python3 docs/reports/gaussian-splitting-small-spacing-sign/check.py` from the repository root. Python 3's standard library is sufficient. The expected final line is `ALL_CHECKS_OK`; the complete output is in `output.txt`.

The source is D. Mikhin and A. Xiourouppa, *A flexible and robust approach to univariate Gaussian splitting using parameterised Gaussian mixtures*, arXiv:2606.01530v1, equations (9)--(11) and (16)--(26). The conjecture following (25)--(26) asks whether the quadratic coefficient of the optimized squared mismatch is negative, for odd component counts. The companion is arXiv:2606.04440v1. Signed symmetric normalized weights are allowed.

The program uses integers and `fractions.Fraction` exclusively. It expands the bilinear forms in equation (10), then computes the source's exact coefficient recurrences for every M from 1 through 8 at sigma = 1/7 and 7/17. It checks both the optimized limiting value against a polynomial Gram projection and the proposed quadratic coefficient against

\[
L^2_2=-\frac{M\Gamma(M+5/2)}{12\pi\sigma(M-1)!}I_{M-1}(a)I_M(a),\qquad
I_n(a)=\int_0^a\frac{t^n}{(1+t)^{n+5/2}}\,dt,\quad a=\frac{1-\sigma^2}{2\sigma^2}.
\]

Multiplication by sqrt(pi) makes all tested coefficients rational. With u = (1+a)^(-1/2), substitution gives the exact polynomial formula

\[
I_n(a)=2\sum_{k=0}^n(-1)^k\binom nk\frac{1-u^{2k+3}}{2k+3}.
\]

For these two widths, u is respectively 1/5 and 7/13. The program also evaluates the source's spatial normal equations directly for M = 1, 2 at s = (h/sigma)^2 equal to 1/1000000 and 1/4000000. All exponentials are enclosed between alternating Taylor sums of degrees 41 and 40. Rational interval elimination encloses the minima, difference quotients and Richardson extrapolations. Each extrapolation differs from the proposed coefficient by less than 1/100000000 of its magnitude. The displayed endpoints are rational outward enclosures; no floating-point computation supplies evidence.

These finite checks verify the source transcription, normalization and algebra in the proposed proof route. They establish neither convergence of the source's formal series nor the universal conjecture. The intended settlement uses a general moving-subspace and Laguerre-polynomial argument, with no finite positive enumeration as its main content. No Lean proof is included here.
