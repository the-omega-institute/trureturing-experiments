# m-resistant pure states for every N and m (Quinta–André–Burchardt–Życzkowski conjecture)

`check.py` (Python with NumPy; exact `Fraction` entries for the reduced states) builds, for N = 3..6 and every 0 ≤ m ≤ N − 2, the pure state
ψ = (2s)^(−1/2) Σ_{S ⊆ [N], |S| = N−m} Σ_{b ∈ {0,1}} ⊗_i |S, b·[i ∈ S]⟩,
with local dimension d = 2·C(N, N−m). It checks two properties:
- every reduced state after losing m parties has a partial transpose with a negative eigenvalue, so it is not fully separable;
- every reduced state after losing m + 1 parties is diagonal in a product basis, so it is fully separable.

`check.log` holds the run (exit 0, `ALL_OK`).
