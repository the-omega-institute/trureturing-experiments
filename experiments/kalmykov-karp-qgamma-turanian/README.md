# Kalmykov–Karp q-gamma Conjecture 1: coefficient clause

Companion programs for the refutation of the coefficient clause of Conjecture 1 in S. I. Kalmykov and D. B. Karp, *Inequalities for series in q-shifted factorials and q-gamma functions*, J. Math. Anal. Appl. 460 (2018) 332–351, arXiv:1702.03788: for doubly positive {g_n}, α > 0, β > 0, µ ≥ max(−β, −1) and α ≤ β + 1, the generalized Turánian of g(µ; x) = Σ g_n x^n/Γ_q(µ + n) is claimed to have positive power series coefficients.

    python3 kk-qgamma-cert.py    # exact rational certificate, g_n = 1, q = (99/100)^10, µ = −7/10, α = β = 7/10, degree 2 (standard library)
    python3 kk-qgamma-scan.py    # mpmath scan of the minimal coefficient over t = α = β, µ, degree ≤ 4, at q = 0.5, 0.9, 0.97
