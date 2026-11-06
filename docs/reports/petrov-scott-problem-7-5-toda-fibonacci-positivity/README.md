# Petrov–Scott Problem 7.5: Toda preservation of Fibonacci positivity

`check.py` is an exact-arithmetic preregistration check for Problem 7.5 of Petrov–Scott, arXiv:2412.21126v2. It reproduces the source's first four Motzkin moment polynomials and the Poisson Toda solution, then checks rational finite determinant, scaling, coordinate-flow, atom, and zero-eigenvector identities used by the proposed refutation family. It uses integers, `fractions.Fraction`, and SymPy exact algebra; no floating-point computation is used as evidence.

The final `ALL_CHECKS_OK` line certifies only these finite algebraic checks. The infinite-support spectral bridge, exponential tilting, and quantifier-level immediate-exit argument remain mathematical work for the post-registration probe.
