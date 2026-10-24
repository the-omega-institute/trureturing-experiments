# Kalmykov–Karp 2017 Conjecture 1: the function clause

Numerical and exact-rational companion to the refutation of the function clause
(Δf(α, β; x) ≥ 0 for x ≥ 0) of Conjecture 1 in S. I. Kalmykov and D. B. Karp,
*Log-concavity and Turán-type inequalities for the generalized hypergeometric
function*, Analysis Math. 43 (2017), arXiv:1611.06724.

Parameters: p = 0, b = (1/2, …, 1/2) of length q, µ = 1, α = β = 1/2.

    python3 kk1f-cert.py    # exact rational certificate at q = 12, x = 5000 (standard library only)
    python3 kk1f-scan.py    # mpmath 80-digit scan of min_x Δf / f(3/2; x)^2 for q = 9..18

The kernel-checked result lives in the main repository at
D5/S3/Analytic/HypergeometricShiftTuran/FunctionResult.lean.
