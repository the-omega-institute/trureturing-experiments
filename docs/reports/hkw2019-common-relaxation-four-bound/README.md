# Common-relaxation four-bound check

Exact checks for Hegde, Keinert and Weber, *A Kaczmarz Algorithm for Solving Tree Based Distributed Systems of Equations*, arXiv:1904.05732v1, Section 3.1 Example 1.

Run with Python 3.12 or newer and SymPy 1.14:

```sh
python3 check.py > output.txt
```

The fidelity check reproduces the source’s two-row averaged example at angle pi/3, including its iteration matrix, eigenvalues and convergence endpoint 8/3.

The candidate uses a rooted star with 2m leaves, weights 1/(2m), root row e0 and leaf rows 4e0 +/- 7ui + 4vi. Every vertex, including the root, uses the same relaxation parameter. For m = 2, 3, 4 the program constructs the original matrices and verifies exact nonsingularity, the affine error equation, a similarity to two scalar blocks and one 2-by-2 block, and the characteristic polynomial. It also verifies the general-m trace, determinant and three quadratic Schur/Jury expressions symbolically. These give the candidate endpoint 1 + 81m/(16(m+1)), whose first three values are 35/8, 307/64 and 101/20.

Each size has 4095 exact rational interior grid checks; just above the proposed endpoint the block determinant exceeds 1. Decimal spectral radii are supplementary diagnostics requested by the preregistration brief; only exact arithmetic supports the mathematical checks. Finite grids do not prove a continuous interval or the arbitrary-m conclusion. The symbolic identities support the separate mathematical argument; no Lean proof or kernel verification is supplied here.

The stored output ends in ALL_CHECKS_OK. The script prints its own SHA-256 digest. This report supports a preregistration draft; completeness of the prior-art screen is separate from these algebra checks.
