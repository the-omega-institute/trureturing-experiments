# Lu–Wenzel complex eigenvalue conjecture: numerical supremum of λ₁ + λ₃

Target: for X ∈ M_n(ℂ) with ‖X‖_F = 1 and T_X(Y) = [X*, [X, Y]] on M_n(ℂ) with the Hilbert–Schmidt inner product, λ₁(T_X) + λ₃(T_X) ≤ 3 (Ge–Li–Lu–Zhou, arXiv:1908.06624; survey arXiv:2402.01085, Conjecture 4.11).

`opt.py n seed` maximizes λ₁ + λ₃ by BFGS from 60 random complex starts (NumPy, SciPy 1.18.1). `run1.log` holds `opt.py n 1` for n = 2, 3, 4, 5:
- n = 2: maximum 2;
- n = 3, 4, 5: maximum 3 to 1e-13, attained by maximizers with different top spectra, all non-normal.

No counterexample was found. The supremum 3 is attained on a set of maximizers with varying spectra.
