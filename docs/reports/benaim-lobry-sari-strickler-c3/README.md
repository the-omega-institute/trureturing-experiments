# Benaïm–Lobry–Sari–Strickler Conjecture 3

`check.py` is an exact symbolic check of the three-season system in arXiv:2407.07553v2, Section 3.2.2. It reproduces the source season order `A1,A3,A2`, verifies the exponential propagator entries, and checks the cyclic reduction

\[
M=\exp(bT)(P^{-1}R)^N,
\qquad R=I+(u-1)E_{11}+vE_{21}.
\]

For every `N` from 3 through 8 it verifies the characteristic polynomial
\(y^N-vy^{N-1}-u\), the proposed positive eigenvector after imposing its root equation, and the exact endpoint expressions used for the large-migration and long-period limits. The script uses SymPy symbols, integers, and rational arithmetic; it uses no floating-point evidence. `check.out` records the run and ends with `ALL_CHECKS_OK`.
