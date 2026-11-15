# Majdoub's localized L6 question

This report checks the proposed moving-soliton obstruction to the small-mass localized L6 estimate in Mohamed Majdoub, *Uniform local smoothing for the derivative NLS and the soliton mass threshold*, arXiv:2611.05330v1, Section 5.5. It supplies preregistration evidence, without a Lean proof or a claim of a completed formal settlement.

Run with Python 3.12 or newer and SymPy 1.14.0:

```sh
python3 check.py
```

`check.out` contains the exact output. The program uses symbolic rational functions with `a = tan(theta)` and `y = exp(x)`. The parameter restrictions are `0 < a < 1`, equivalently `0 < theta < pi/4`. It verifies the complex travelling-wave PDE, the amplitude and phase equations, mass `8*theta`, scaling covariance, and the transit exponent. Its independent source-fidelity check reproduces the source's integral `W^2 = 512/7`. The rectangle inclusion identity supports the explicit positive linear lower bound for the localized L6 integral; the inequalities and change of variables are explained in the preregistration draft. Symbolic computations check these identities; they do not certify the full analytic argument in Lean.

For additional quadrature diagnostics, install mpmath 1.3.0 and run:

```sh
python3 check.py --numerical
```

`numerical.out` contains that output. Floating point values are diagnostics only and are excluded from the exact acceptance checks. The computation uses the exact polynomial primitive of `sech(x/99)^36` and integrates over the soliton profile. The displayed `I6/lambda` approaches the positive predicted constant at two fixed masses, one below 1 and the other equal to pi. The last line `ALL_CHECKS_OK` reports the exact checks, and makes no numerical-certification claim.

The check program is original mathematical verification code. SymPy and mpmath have BSD-style licenses; no external program or Lean proof was copied. The soliton formula and the two thresholds are attributed to Majdoub and to Harrop-Griffiths–Killip–Ntekoume–Vişan, DOI 10.4171/JEMS/1490, Section 1.3.
