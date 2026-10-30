# Sedlák–Reitzner–Chiribella–Ziman maximal tester orthogonality

`check.py` uses only `fractions.Fraction` arithmetic. It checks the source's
orthogonal qubit testers and exhaustively checks, for small finite dimensions,
outcome counts, and rational diagonal tester effects, the exact diagonal
covering identity

\[
t^* = \sum_x \max(\rho_x,\sigma_x),\qquad R_t=1-1/t^*.
\]

The final assertion checks `R_t = 1/2` exactly when the two normalization
supports are disjoint, equivalently when every diagonal effect product is zero.
