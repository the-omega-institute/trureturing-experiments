# Petrov–Scott, Conjecture 6.3: the alternative q-Charlier specialization is Fibonacci positive

Check for the open-problem preregistration proving Conjecture 6.3 of L. Petrov and J. Scott, "Random Fibonacci words via clone Schur functions", Forum of Mathematics, Sigma 14 (2026), e15 (arXiv:2412.21126v2).

Run: `python3 check.py` (Python 3 with sympy; exact rationals).

For the alternative q-Charlier specialization x_k, y_k, its normalization c_k, t_k (Proposition 4.4 of the source), the tridiagonal determinants A_ℓ, B_k of (2.10) and the clone Schur functions s_w of Definition 2.5, the program checks:
- symbolically for k = 1, …, 8: x_k = c_k(1 + t_{k−1}) with t_0 = 0, and y_k = c_k c_{k+1} t_k;
- at ρ = q = 1/4, by exact determinants: s_w > 0 for all 609 words w over {1, 2} with |w| ≤ 12;
- at ρ = q = 1/4: 0 < t_k ≤ 1/4 and t_{k+1} ≤ (125/256) t_k for k ≤ 58, and E_k(m) ≥ (143/768) t_{m+1} for k, m ≤ 20, where E_k(m) = t_{m+1} − (1 + t_m − t_{m+1}) t_{m+2} D_{k−1}(m+3) and D_n(m) = 1 + Σ_{r=1}^{n} ∏_{j<r} t_{m+j};
- for k, ℓ, m < 8: B_k(m) = E_k(m) ∏_{j=m+1}^{m+k+2} c_j and A_ℓ(m) = D_ℓ(m) ∏_{j=m+1}^{m+ℓ} c_j.

`run.log` holds the output of the run cited by the preregistration.
