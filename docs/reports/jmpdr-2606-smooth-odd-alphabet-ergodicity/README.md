# `jmpdr2606p2` exact checks

`check.py` implements the type-1 substitution from Jamet–Marcovici–Poirier–de la Rue, specialized to `α=1` and `β=2q+1`. It decodes `A,B,C,D` as `1^1, β^1, 1^β, β^β`, exhaustively checks the run-length identity on every alternating word of lengths 1 through 8, and verifies the incidence matrix, characteristic polynomial, normalized eigenvectors, and distinct decoded frequency limits with exact SymPy arithmetic.

The integration check covers `q=1,…,8` (8,160 words); the requested refutation family is the strict-gap range `q≥3`, where `q>λ=(1+√(4q+1))/2`. The boundary `q=1` is the source `{1,3}` theorem, while `q=2` is the spectral-tie boundary and is left open here. The script uses no floating-point arithmetic and ends with `ALL_CHECKS_OK` on success.
