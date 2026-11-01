# Finch third primitive determinant

Exact checks for Steven Finch, *Small Matrices with Large Inverses: Unimodular 4 by 4 Cases*, arXiv:2607.07688v4, Section 4 (Closing Words).

Run with Python 3 and SymPy:

```sh
python3 check.py
```

`output.txt` contains the output. The program reproduces the source's explicit extremal determinant and cube ceiling, and its first three proposed primitive levels for k = 3 through 7. It additionally checks k = 8, 9. For each k it visits every ordered pair of rows; the exact upper bound k times the sum of positive cofactors excludes all omitted matrices from the range at or above the proposed third level. It computes all nine minors for primitivity. Witness bounds, determinants, gcds and ordering are checked through k = 1000.

The symbolic checks generate 28 decoded row pairs and twelve cofactor types up to interchanging the two positive coordinates, then check polynomial identities and loss inequalities using exact integer coefficients. These computations support the preregistered general argument. They are not a Lean proof or a claim of settlement for every k.

Dependencies: Python 3 standard library and SymPy 1.14.0. No floating point is used as evidence. Successful output ends with `ALL_CHECKS_OK`.
