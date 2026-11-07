# Magnetic billiard rotational rigidity

Source: Michael (Misha) Bialy, *Magnetic billiards: duality and magnetic Radon property*, arXiv:2610.33782v1, Section 13, question 5.

Run `python3 check.py` with Python 3 and SymPy 1.14.0. The stored standard output is `check.out`. Assertions use exact integers, rationals and symbolic algebra; no floating point enters the checks.

The program reproduces the source's Section 4 endpoint identities and circular factorization. It verifies the proposed radial velocity-jump identity symbolically and tests all nonzero rotation steps at orders 2, 3, 4 and 6, for two table radii and three larger Larmor radii each. The configurations satisfy the two equal-radius constraints and specular reflection. An exact ellipse example tests that radiality of the normal is a substantive restriction.

Solving the source's Section 3 linear system gives `t'_beta = rho_beta sin(s)/(R sin(mu))`. The source's subsequent display has the opposite sign and omits `1/R`; its displayed mixed derivative also retains an extra factor `R`. The corrected mixed derivative is `-rho_alpha rho_beta sin(t) sin(s)/sin(mu)`. Its asserted negative sign, and therefore the qualitative twist conclusion, agree with the corrected calculation.

These checks support the local algebraic mechanism. They do not prove global collision-map existence, the Birkhoff graph theorem, the all-periodic graph uniqueness theorem, or the full conjecture.
