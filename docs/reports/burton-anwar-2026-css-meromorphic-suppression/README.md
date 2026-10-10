# Burton–Anwar arXiv:2605.06251 Conjecture 4.8

Numerical check for the settlement preregistered in https://github.com/the-omega-institute/trureturing/issues/14813. The kernel-checked result and its dossier `Problems/burton-anwar-2026-css-meromorphic-suppression.md` are in the main repository; this directory holds only the experimental program.

Run from this directory: `python3 check.py (exact arithmetic, standard library)`.

Files:
- `check.py`: sha256 `3a59a3e7d8503363420a9c0ed3738e8d611a06de9e74effc0caa6901953a684b`

`check-source-example.py` (SymPy) computes the local degrees of the source's [[15,1,3]] decoder f(z) = (z^15 + 15z^7)/(15z^8 + 1) at 0, ∞, 1 and −1; `check-source-example.log` holds its run (exit 0; degrees (7,7,3,3)).
