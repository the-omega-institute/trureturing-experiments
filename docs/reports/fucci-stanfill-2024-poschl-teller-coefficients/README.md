# Fucci–Stanfill: a vanishing coefficient g_{q,1}(p/q, β) = 0

`check.py` checks a family refutation of the conjecture in Remark B.3 of arXiv:2411.17860v1 (Ann. Henri Poincaré 27 (2026) 2073–2116): "all of the g_{m,j}(p/q,β) are nonzero".

- It rebuilds P₁(y,z) from the source's general definitions: E_k and C_m with Bernoulli numbers, and G_k(x,y) = binom(x−y,k) B_k^{(x−y+1)}(x) with the generalized Bernoulli polynomials of DLMF 24.16.1. It compares the result with the printed formula (B.21) and checks that the coefficient of sin α/[…] is symmetric under y ↦ 1−y.
- For every coprime 0 < p < q ≤ 8, with Ω₀ = C, Ω₁ = C[P₁((1+ν)/2) − P₁((1−ν)/2)] and arbitrary symbolic higher Ω_k, it extracts S_q, the coefficient of x^{q+p} in log(1 + Σ Ω_k x^{p+kq}). It checks that g_{q,1} = 0, that S_q does not depend on the higher Ω_k, and that g_{q,0} = Cν(ν²−1)/3 for p ≥ 2.

Run: `python3 check.py` with sympy (exit 0, final line `ALL_OK`).
