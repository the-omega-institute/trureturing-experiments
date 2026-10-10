# Allouche–Shallit Conjecture 12: sixfold evil/odious representation counts

`check.py` checks the sharp clauses of Conjecture 12 of arXiv:2112.13627v3 and the prefix-certificate proof route.

- Direct sixfold convolutions for n < 8192: r₆ is strictly increasing from 37 and s₆ from 5; r₆(36) > r₆(37).
- The first-difference identities 64Δr₆(n) = C(n+4,4) + Σ_{j=1}^{5} C(6,j) h_j(n) + c(n) − c(n−1) and the analogue for s₆ with signs (−1)^j, where h_j and c are the coefficients of T^j/(1−z)^{5−j} and T^6, T(z) = Σ(−1)^{t(n)} zⁿ.
- The two-digit state matrices of the dyadic recursions, their row-sum norms, and the exhaustive length-5 products of the T^6 system.
- The 12-bit prefix certificate: p⁴ − 24·B_p > 0 for 2048 ≤ p ≤ 4095.

Run: `python3 check.py` (exit 0, final line `ALL_OK`).
