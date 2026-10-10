# Hirche–Guan–Tomamichel: order-three information combining with quantum side information

`check.py` checks numerically the order-three equality of arXiv:2305.02589v1 (the remark after Proposition V.5 and the α = 3 equality clause of Conjecture V.5): for independent binary classical–quantum inputs ρ_i = ½|0⟩⟨0|⊗σ_i0 + ½|1⟩⟨1|⊗σ_i1 and their XOR combination τ,

H̃₃↓(X₁+X₂|B₁B₂)_τ = −½ log[(4K₁K₂ − K₁ − K₂ + 1)/3],  K_i = exp(−2 H̃₃↓(X_i|B_i)),

where H̃₃↓(A|B)_ρ = −½ log Tr[(ρ_B^{−1/3} ρ_AB ρ_B^{−1/3})³] with inverse powers on the support of ρ_B. It also checks that this value equals both branches of the conjectured BSC-PSC expression at α = 3 (with h₃(p) = −½ log(p³+(1−p)³), its inverse on [0, ½], and binary convolution).

- 400 random trials, output dimensions 1–4, random ranks (including rank-deficient σ's), noncommuting σ's.

Run: `python3 check.py 7` with numpy (exit 0; final line `trials=400 max_abs_err=1.343e-14`).
