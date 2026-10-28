# Alhasan–Fawwaz density-independent bound: exact checks

Source: S. A. Abed, I. Nikolaeva and A. Novikov, *How characterisations of traces and central elements fail: sharp violation constants via Michelson contrast and noncentrality*, arXiv:2610.34388v1, Proposition 5.44(d) and Section 11.2(5).

Run with Python 3 and SymPy 1.14.0:

```sh
python3 check.py
```

The stored output is `output.txt`. All evidence uses exact integer and rational complex arithmetic. The program reproduces the source's three-by-three example, including its displayed principal minors and weighted quotient. It checks a rational-coordinate triangular sharpness construction in dimensions 1 through 5, with positive definite modulus and density, quotients exceeding `n - 1/100`, and a singular extension with the claimed rank support. It also checks every accretive input in the explicitly enumerated two-by-two grid: 320 candidates with positive definite modulus, of which 224 are accretive.

For the construction, set `c_j = (1-j²t²)/(1+j²t²)`, `lambda_j = (1+ijt)/(1-ijt)` and let `T` have diagonal `lambda_j` and upper entries `2c_k`. With `W = diag(1/c_j)` one has `(TW + WT*)/2 = ee*`. Construct a unit upper triangular eigenvector matrix `S` with `TS = SD`, and set `R = S⁻¹`, `P = RWR*`, `X = DP`. Then `|X| = P` and `Re X = y y*`, where `y = R e`. The checker verifies these identities exactly and supplies a strictly positive density whose weighted quotient exceeds the displayed threshold.

These finite checks support the preregistration's proposed proof route. They do not establish the universally quantified assertion or constitute a Lean proof.
