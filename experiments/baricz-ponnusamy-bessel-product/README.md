# Baricz–Ponnusamy conjecture: log-convexity of ν ↦ I_ν(u)K_ν(u)

Companion programs for the refutation of the conjecture at the end of Á. Baricz and S. Ponnusamy, *On Turán type inequalities for modified Bessel functions*, Proc. Amer. Math. Soc. 141 (2013) 523–532, arXiv:1010.3346: ν ↦ I_ν(u)K_ν(u) is log-convex on (−1, ∞) for every u > 0.

    python3 bp-cert.py    # exact rational certificate at u = 1, orders −1/2, 1/2, 3/2 (standard library)
    python3 bp-scan.py    # mpmath check of the closed forms and of midpoint log-convexity gaps
