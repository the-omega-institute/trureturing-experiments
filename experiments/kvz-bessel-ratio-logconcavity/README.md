# Karp–Vishnyakova–Zhang Conjecture 1: log-concavity clause

Companion programs for the refutation of the log-concavity clause of Conjecture 1 in D. Karp, A. Vishnyakova and Y. Zhang, *Unimodality preservation by ratios of functional series and integral transforms*, Results Math. (2025), arXiv:2408.01755: the ratio x ↦ I_{ν1}(a1 x)/I_{ν2}(a2 x) is claimed log-concave on (0, ∞) for ν1 ≥ ν2 > 0, 0 < a1 ≤ a2.

    python3 kvz-cert.py                      # exact rational certificate, ν1 = ν2 = 1/10, a1 = 1, a2 = 2, t = 3/2, 9/2, 15/2 (standard library)
    python3 kvz-scan-second-derivative.py    # mpmath scan of d²/dx² log ratio over (ν1, ν2, a2)
    python3 kvz-scan-midpoint.py             # mpmath scan of midpoint log-convexity gaps
