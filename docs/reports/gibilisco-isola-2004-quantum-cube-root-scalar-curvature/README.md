# Gibilisco–Isola: monotonicity of scalar curvature for the quantum cube-root (p = 3) α-geometry

`check.py` supports a proof of the p = 3 slice of the second clause of Conjecture 4.2 of math-ph/0407007v2 (J. Math. Phys. 46 (2005) 023501): on the positive definite density matrices with the pull-back of ρ ↦ 3ρ^(1/3) under the real Hilbert–Schmidt inner product, Scal_3 is strictly Schur-increasing for every n ≥ 2. With y_i = λ_i^(1/3) and s_k = Σ y_i^k, the formula checked is Scal_3 = N/(9 s_4²), N = [(4n²−2)s_1² − 2n s_2]s_4 − 8n s_1 s_5 + 8 s_6.

- (1) Intrinsic scalar curvature of the pull-back metric in coordinates (trace-orthogonal Hermitian basis, Daleckii–Krein derivative of the embedding, Christoffel symbols and Ricci contraction by nested central differences at 60 digits), at one point for n = 2 and one for n = 3, compared with N/(9 s_4²). Control: the same code gives (n²−1)(n²−2)/4 for p = 2 (the source's value).
- (2) The Gauss-equation value (tr S)² − tr(S²) of the hypersurface Tr Y³ = 1, with the shape operator built as an explicit n² × n² matrix at random unitary conjugates, three points each for n = 2..5, compared with N/s_4².
- (3) Symbolic algebra: the chain-rule identity for the two-coordinate derivative (n = 2, 3, 4); divisibility defining Q; the moment-deficit identity with symbolic n; the closed form of L; positivity of L for n ≥ 4 (all 126 coefficients of L(n = k+4) positive), for n = 2 (11 positive coefficients) and for n = 3 (41 positive coefficients after subtracting y²z⁴(2x² − z²)²).
- (4) Strict increase of N/(9 s_4²) under 19,999 random nontrivial T-transforms (n = 2..7).

Run: `python3 check.py` with sympy 1.14 and mpmath 1.3 (exit 0; final line `bad= 0`; readings `n=2 intrinsic Scal_3 = 1.901978716728070992646611`, `n=3 intrinsic Scal_3 = 16.77929912505800523121669`, each equal to N/(9 s_4²) to the printed digits).
