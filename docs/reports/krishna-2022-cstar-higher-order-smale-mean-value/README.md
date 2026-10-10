# Krishna arXiv:2206.08154 Higher Order C*-algebraic Smale Mean Value Conjecture: counterexample in ℂ × ℂ

`check.py` (SymPy, exact roots and critical points) works in A = ℂ × ℂ with the supremum norm, at n = 3, z = 0 and P = (x³ + 21x² + x, y³ − 3y). It checks:
- ‖P′(0)‖ = 3 and ‖P″(0)‖ = 42;
- for each of the four critical points w, the k = 2 quantity (‖P″(z)‖/2!)·‖P(z) − P(w)‖/‖P′(z)‖² exceeds 4, with a minimum of 14/3.

For general t > 18 the second coordinate gives |P₂(0) − P₂(±1)| = 2, so the k = 2 quantity is at least 2t/9 > 4. `check.log` holds the run (exit 0).
