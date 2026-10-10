# Krishna arXiv:2206.06653 Conjectures 2.1, 2.3, 2.4

Numerical check for the settlement preregistered in https://github.com/the-omega-institute/trureturing/issues/14845. The kernel-checked results and their dossier are in the main repository; this directory holds only the experimental program.

Run from this directory: `python3 check.py` (SymPy).

It verifies, exactly and with symbolic `t, c > 0`, the factorization of the ordered derivative for the M₂(ℂ) witness family and the six defect matrices of the three conjectures.

Files:
- `check.py`: sha256 `702b68b36e2ae270381994118ce00660b5310f00ad5c514705b95334a7b0f4dd`

`scalar-kushel-tyaglov.py` (SymPy, exact) evaluates the scalar inequality that the source quotes as the Kushel–Tyaglov theorem (Section 1), at d = 3 and a = (1+i, 1−i, 1). It checks the factorization of P′, gets lhs = 32/9 and rhs = 28/9, and so finds the inequality false by 4/9. `scalar-kushel-tyaglov.log` holds its run (exit 0).
