# Diaconis–Miclo Conjecture 1.3

`check.py` uses rational interval arithmetic and SymPy exact polynomial arithmetic. It verifies the five cyclotomic character parameters, the literal 4×4 Fourier-block determinant, strict Schur certificates at `a = 1/5` and `a = 3/10` inside radius `239/500`, and an exact sign change at `a = 1/4` in the middle block. The sign change places a real eigenvalue below `-239/500`; the two endpoint certificates place every nontrivial one-coordinate block inside that radius. No floating-point value is used as evidence.

The program is an independent arithmetic check, not a Lean proof and not a complete formalization of the all-dimensional product-family argument.
