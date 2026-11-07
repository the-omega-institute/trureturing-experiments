# Bounded-class quadratic sharpness on se(3)

The check concerns the unnumbered bounded-class sharpness question after the comparison of Proposition 3.4 and Theorem 4.2, repeated in Section 6 of Sooraj K.C and Vivek Mishra, *Operator-norm bounds and a quadratic lower-growth example for the special Euclidean algebra se(3)*, [arXiv:2605.31076v1](https://arxiv.org/abs/2605.31076v1).

Run with Python 3 and SymPy:

```sh
python3 docs/reports/kcm2605se3-quadratic-sharpness/check.py
```

The recorded output uses Python 3.14.8 and SymPy 1.14.0. Every calculation uses exact symbolic expressions, integers or rational numbers; no floating point supplies evidence.

The fidelity checks reproduce the source's distinction between the six-vector norm and the matrix Frobenius norm, and its two Rodrigues coefficient limits. A fresh matrix exponential computation on the rotation/translation family gives its translation entry and the exact second-order Taylor coefficient of the proposed cosine objective. Ambient differentiation gives gradient norm at most one and Hessian operator norm at most one. Exact integration checks the coefficients in the proposed derivative estimates, and rational radius samples check the resulting inequalities.

The proposed fixed objective is F(G) = cos(G_24), with one-based matrix indices. The analytic outline gives R²/16 ≤ L_F(R) ≤ 25R²/12 for R ≥ 1. The program checks symbolic identities and finite rational samples; the global Duhamel estimates and passage from a directional second derivative to the Lipschitz supremum remain analytic proof obligations. This experiment is neither a Lean proof nor an independently certified settlement.

Files: `check.py` contains the check; `output.txt` contains its complete successful output, ending with ALL_CHECKS_OK.
