# Cotfas (arXiv:2610.05623): positivity of the mixed discrete-variable Gaussian state

Source: after Eq. (22), "As concern ρ̃_σ, we do not know if (11) is sufficient in order to have ρ̃_σ ≥ 0. Based on our numerical investigations, we conjecture that this is true."

Route (proof sketch from the search seat, to be kernel-checked):
- **Reduction to h = d.** Rescaling the arguments of W gives M_h = (d/h)·M_d.
- **Kernel identity.** For h = d, the matrix satisfies (M_d)_{ij} = 2 Σ_{r,l∈ℤ} K(i + dr, j + dl), where K(x, y) = g(x)·conj(g(y))·e^{−λ(x−y)²}, g is a complex Gaussian and λ = π(det σ − 1)/(2da) ≥ 0. The identity follows from Poisson summation on ½ℤ.
- **Positivity.** e^{−λu²} is the characteristic function of a centred Gaussian (a point mass when λ = 0). This writes M_d as a positive mixture of rank-one matrices.
- **Positive trace.** A nonvanishing theta value gives tr M_d > 0, so normalization is defined.

## Programs

Both need numpy (`python3 -m pip install numpy`).
- **`check.py`.** Run `python3 check.py`; it prints `RESULT PASS` when every check holds. It covers 600 covariances with det σ ≥ 1 over d ∈ {1, 3, 5, 7, 9} and h ∈ {0.5, 1, 3}. Each matrix must be Hermitian and positive semidefinite, with positive trace. It also reports the 36 sampled covariances with det σ < 1.
- **`identity.py`.** Run `python3 identity.py`. It prints the largest entrywise error of the kernel identity for d ∈ {3, 5, 7}.

Readings (`run.log`, exit 0):
- **`check.py`.** 600 cases, 0 failures; the most negative relative eigenvalue is −7.7·10⁻¹⁶ (rounding). All 36 cases with det σ < 1 are not positive semidefinite.
- **`identity.py`.** Largest entrywise error 1.7·10⁻¹⁵.
