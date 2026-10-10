# Bunth–Pitrik–Titkos–Virosztek: the quantum Wasserstein complexity of channels is not convex

`check.py` addresses the open question of arXiv:2402.13150v4, §6.1 ("whether the Wasserstein complexity … is convex or not is an open question"), with C_W(Φ) = max_ρ d_A(ρ, Φ(ρ)) and the self-cost-corrected quantum Wasserstein divergence d_A of De Palma–Trevisan couplings.

For a qubit with A = (X, Y, Z), t = 1/4, the unitary channels Φ_± = Ad(√(1−t) I ± i√t Z) and their midpoint Ψ = (1−t) id + t Ad(Z), it solves the coupling semidefinite programs (cvxpy) and checks:
- the self cost D²(σ,σ) = 8 − 4 (tr √σ)² on random states;
- d²(ρ, Φ_±(ρ)) ≤ 8t = 2 on 300 random states (the bound is proved for every state by an explicit coupling);
- d²(|+⟩⟨+|, Ψ(|+⟩⟨+|)) = 1 + √3 > 2.

Hence C_W(Ψ) ≥ √(1+√3) > √2 ≥ (C_W(Φ_+) + C_W(Φ_−))/2.

Run: `python3 check.py` with numpy and cvxpy (exit 0; output lines `self-cost formula ok: True`, `max d^2 over sampled rho for Phi_+-: 1.654607 (bound 8t = 2)`, `d^2(+, Psi(+)) = 2.732051  1+sqrt3 = 2.732051`, `bad= 0`).
