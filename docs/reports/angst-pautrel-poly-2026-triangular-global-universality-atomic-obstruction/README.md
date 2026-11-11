# Triangular signals: an atomic obstruction to global expected zero counts

Source: Jürgen Angst, Thibault Pautrel and Guillaume Poly, *Global universality of the expected number of zeros of non-analytic random signals*, [arXiv:2609.01007v1](https://arxiv.org/abs/2609.01007v1), Remark 5(2). The source counts the entire zero set and permits iid centered unit-variance coefficients with a finite third moment.

For the source's continuous periodic triangle, the repeated coefficient block `(1,-1,-1,1)` has zero zeroth and first index moments. At degree `n=4m`, its signal is zero throughout `[0,π/n]`, while its value at `π/2` is `2m`. An iid Rademacher sequence realizes this prefix with probability `2^(-4m)>0`. Thus the full zero-set cardinality, and its extended expectation, are infinite on these degrees. More generally the same mechanism applies to any admissible law with positive masses at `c` and `-c`, for `c>0`; the event probability is `(p_+ p_-)^(2m)`. This addresses the unrestricted coefficient class. It gives no settlement for purely nonatomic laws, including Gaussian coefficients.

The check reproduces the triangle's interpolation of cosine on `πℤ` and the cosine normalization in Remark 5(1). Exact symbolic integration gives the triangle's proposed constant `4/π`. It verifies the symbolic block cancellations and checks degrees `4m` for `m=1,...,100`, three amplitudes and five rational interval locations. An independent enumeration counts every sign vector whose first cell is identically zero at degrees 4, 8 and 12. Finite checks support the displayed algebra; the infinite-family conclusion requires the mathematical argument, not extrapolation from enumeration. No Lean proof is included.

Run with Python 3.12 or newer and SymPy 1.14:

```sh
python3 docs/reports/angst-pautrel-poly-2026-triangular-global-universality-atomic-obstruction/check.py
```

All evidence uses integers, fractions and exact symbolic expressions. The stored output is `check.out`; its final line is `ALL_CHECKS_OK`.
