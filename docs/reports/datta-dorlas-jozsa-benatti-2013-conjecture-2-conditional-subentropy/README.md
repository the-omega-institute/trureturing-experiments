# Datta–Dorlas–Jozsa–Benatti (arXiv:1310.1312) Conjecture 2: conditional subentropy below the min-entropy

Source: N. Datta, T. Dorlas, R. Jozsa and F. Benatti, "Properties of subentropy", arXiv:1310.1312v1 (J. Math. Phys. 55, 062203 (2014)), Equations (1)–(4), (13)–(17) and Conjecture 2.

Conjecture 2: for a c-q state ρ_XB, Q(X|B) = Q(ρ_XB) − Q(ρ_B) ≤ H_min(X|B) = −ln p_guess(X|B).

The proposed proof gives the stronger bound Q(X|B) ≤ 1 − p_guess(X|B) ≤ −ln p_guess(X|B).

## Program

`check.py` needs numpy, mpmath and cvxpy with the CLARABEL solver (for example `python3 -m venv v && v/bin/pip install numpy mpmath cvxpy`). Q is computed from Equation (1) in 60-digit arithmetic. The program prints `RESULT PASS` when every check holds:
1. **Equation (4).** It agrees with Equation (1) for three spectra with n = 2 and 3 (error below 1e−31).
2. **120 random c-q states.** They have dim B = 2–4 and 2–4 labels, mixed and singular blocks, and a third of them nearly distinguishable.
   - p_guess is bracketed by a feasible POVM and a dual-feasible σ ≥ A_x, so P_up = tr σ ≥ p_guess. The largest primal–dual gap is 2.3e−8.
   - Every state has 0 ≤ Q(X|B) ≤ 1 − P_up, which implies Q(X|B) ≤ 1 − p_guess. The largest value of Q(X|B) − (1 − P_up) is −0.0019.
3. **Lemma 1.** Q(r) ≤ 1 − max_j r_j holds for 400 random probability vectors (largest difference −6.2e−9, near one-hot vectors). The one-sided derivative f′(0+) of f(t) = Q(1 − t, t/(m − 1), …) is 1 for m = 2, 3, 5 and 8.
4. **Lemma 2 (statistical).** Q(ΔR) − Q(R) = D_KL(q_R ∥ q_ΔR) for the pinching Δ by three projections in dimension 4. The check uses Monte Carlo over Haar-random unit vectors (400000 samples, tolerance 5 standard errors).

Reading: `run.log` (exit 0, `RESULT PASS`, 4.6 s).
