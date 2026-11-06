# Sample-covariance logarithmic-energy check

`check.py` uses only exact integer, SymPy rational, and `fractions.Fraction`
arithmetic. It reproduces the square Laguerre formula from Chafaï–Dadoun–Youssef at
`n = 1` and its exact Marchenko–Pastur limit `3/2`. It then enumerates the finite
Bernoulli support for `n = 2` and `n = 3`. For each pair of independent matrices it
factors the Gram characteristic polynomials exactly and computes the small-noise
determinant order

`k(A,C) = 2 z(A) z(C) + Σλ m_A(λ)m_C(λ)`

over common positive Gram eigenvalues. The resulting pairwise-energy coefficients
`L₂` and `L₃` are exact fractions and satisfy `L₃ > L₂`, as well as the bounds used
in the preregistration. This is an arithmetic check of the proposed mechanism; it
does not replace the analytic Gaussian-polynomial integrability and limiting lemmas.
