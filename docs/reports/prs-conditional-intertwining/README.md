# PRS conditional intertwining check

`check.py` is an exact-arithmetic finite reproduction for Adrian Padellaro, Sanjaye Ramgoolam, and Ryo Suzuki, arXiv:2610.10745v1, Section 4.3, Conjecture 2. It constructs the finite permutation centraliser algebra, the Young-height cutoff as the kernel orthogonal complement, the derangement orbit span, their intersection, the coefficientwise orthogonal projection, the permutation operator at ξ = 1/N, and the actual traceless contraction.

The check first reproduces the published traceful identity qH = Δq at (m,n,N) = (2,2,2). It then checks the conditional traceless identity F(PHu) = Δ′F(u) on every basis vector of U for (m,n,N) equal to (1,1,2), (2,2,2), (2,2,3), and (3,1,2). All coefficients use integers or exact `sympy.Rational` values; the script uses no floating point arithmetic. The output is stored in `check.out` and ends with `ALL_CHECKS_OK`.
