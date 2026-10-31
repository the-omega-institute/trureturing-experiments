# Namkung–Kwon singular-boundary residual check

This report supports a preregistration draft concerning Conjecture 2 of Min
Namkung and Younghun Kwon, “Generalized sequential state discrimination for
multiparty QKD and its optical implementation”, Scientific Reports 10, 8247
(2020), DOI [10.1038/s41598-020-63719-9](https://doi.org/10.1038/s41598-020-63719-9).
The target is the intended determinant-zero boundary, allowing zero individual
success probabilities at the second stage.

Run with Python 3 and SymPy:

```sh
python3 docs/reports/namkung-kwon-singular-boundary-residual/check.py
```

The recorded run uses Python 3.14.8 and SymPy 1.14.0. All calculations use exact
integers, rationals and algebraic expressions. A successful run exits zero and
ends with `ALL_CHECKS_OK`; its complete output is in [output.txt](output.txt).

The source-fidelity check reproduces the binary residual Gram of all ones at
equal boundary success and checks its symbolic determinant obstruction, with an
additional exact grid check. The family checks cover n = 4 through 8, symbolic
positive-definite Gram lifts, the dependency supported on the first three
labels, and the identities proving the residual region is
{0}^3 × [0, 1/3]^(n−3). They check 62 box vertices and three rational lift
parameters at each dimension. A separate n = 4, p = 1/2 calculation constructs
exact normalized inputs and Kraus matrices and verifies completeness, zero error,
outcome-independent outputs and a nonzero second-stage unambiguous measurement.

Finite computations do not prove the assertion for every n and p. The draft
supplies the general ordinary-mathematical argument. These checks are neither a
Lean proof nor a claim that literature screening is exhaustive. The source is
quoted and cited in the draft; this directory contains an independently written
check program and its output, with no copied source code or article figures.
