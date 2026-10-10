# Gibilisco–Isola: monotonicity of scalar curvature for the logarithmic (p = ∞) α-geometry

`check.py` supports a proof of the p = +∞ clause of Conjecture 4.1 of math-ph/0407007v2 (J. Math. Phys. 46 (2005) 023501): on the positive simplex with the pull-back of ρ ↦ log ρ, the scalar curvature Scal_∞ is strictly Schur-increasing for n > 2.

- It computes the intrinsic scalar curvature from the metric in coordinates (Christoffel symbols, Ricci contraction; sympy, exact rational points) for n = 3, 4, 5 and compares it with R = 2(e₃ − 4e₄)/s₂². As a control, the same code gives (n−1)(n−2)/4 for p = 2 (the source's sphere value).
- It checks symbolically the identity K = Σ z_i²(1−2z_i)² + u(3−4u)B + u²(1−u) + v(4B − 4u² + 6u − 2), where K = (A − 4E₂)T + 4N is the numerator of dR/dv along two-coordinate mixing, for 1–4 remaining coordinates.
- It checks strict increase of R under 20,000 random T-transforms (n = 3..8).

Run: `python3 check.py` with sympy (exit 0; final line `bad= 0`; readings `R 150/361`, `202/225`, `3712/2025` equal to the prediction at the three sample points).
