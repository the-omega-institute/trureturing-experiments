# Wong–Bravyi–Gosset–Liu, Section 5.6: Conjecture 2 fails for a phase-shifted 4-cycle

Check for the open-problem preregistration refuting the phase-shifted extension of Conjecture 2 of B. Wong, S. Bravyi, D. Gosset and Y. Liu, "Lee-Yang tensors and Hamiltonian complexity", arXiv:2602.03605v1 (2026).

Run: `python3 check.py` (Python 3 with numpy and sympy).

For the 4-cycle with edge phases (0, 0, 0, π) and s = 3/5, H_s(θ) = −Σ_e h_e(s e^{iθ_e}) with h(z) = |00⟩⟨00| + z|11⟩⟨00| + z*|00⟩⟨11| + |z|²|11⟩⟨11|. The program checks, exactly with sympy where stated:
- H is Hermitian and preserves Hamming-weight parity (exact);
- −α₊, −β₊ are roots of the even-block characteristic polynomial and −ω₊ of the odd one, with t = s², α₊ = (5 + t + √((1+t)(9+t)))/2, β₊ = (1 + 5t + √((1+t)(1+9t)))/2, ω₊ = 1 + t + √(1+t²) (exact);
- β₊ − ω₊ = (1 + √901 − √706)/25 > 0 (exact), so the second smallest eigenvalue −β₊ (even, simple) lies below the odd minimum −ω₊;
- the unphased 4-cycle at s = 0.6 satisfies Conjecture 2, and the π-flux 4-cycle violates it for every s ≥ 0.6 on a grid of step 0.05 (and already at 0.55).

`run.log` holds the output of the run cited by the preregistration.
