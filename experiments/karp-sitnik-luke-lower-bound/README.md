# Karp–Sitnik Conjecture 1: lower bound for (q+1)Fq(σ, a; b; −x)

Companion programs for the refutation of Conjecture 1 in D. Karp and S. M. Sitnik, *Inequalities and monotonicity of ratios for generalized hypergeometric function*, J. Approx. Theory 161 (2009) 337–352, arXiv:math/0703084: the lower bound (1 + x∏aᵢ/bᵢ)^(−σ) < (q+1)Fq(σ, a; b; −x) is claimed for all σ > 0 under Σ(bᵢ − aᵢ) > 0.

    python3 ks-cert.py    # exact rational certificate, σ = 1, a = (1, 1), b = (1/10, 2), x = 1/10 (standard library)
    python3 ks-scan.py    # mpmath scan over q = 2 parameter grids with mixed signs of bᵢ − aᵢ
