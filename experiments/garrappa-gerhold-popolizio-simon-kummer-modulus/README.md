# Garrappa–Gerhold–Popolizio–Simon Remark 6: does |M(α, α−β, z)| ≥ M(α, α−β, ℜz) hold for all α, β > 0?

Companion program for the negative answer to the question in Remark 6 of R. Garrappa, S. Gerhold, M. Popolizio and T. Simon, *On some inequalities for the two-parameter Mittag-Leffler function in the complex plane*, J. Math. Anal. Appl. 551(1) (2025) 129588, arXiv:2410.11852. M is Kummer's function ₁F₁.

    python3 ggps-check.py    # mpmath (50 digits): three witnesses where the inequality fails

- A: α = 1/2, β = 2. M(1/2, −3/2, z) = e^z (1 − 4z/3 + 4z²/3) vanishes at z = 1/2 + i√2/2 while M(1/2, −3/2, 1/2) = (2/3)e^{1/2} > 0.
- C: α = 3/2, β = 3. M(3/2, −3/2, z) = e^z (1 − 2z + 4z² + 8z³/3) has a non-real zero with positive value at its real part.
- B: α = 2, β = 3/2 (second parameter 1/2 > 0). At z = −10 + 5i, |M| = 0.010077 < M(2, 1/2, −10) = 0.016389.

Requires `mpmath`.
