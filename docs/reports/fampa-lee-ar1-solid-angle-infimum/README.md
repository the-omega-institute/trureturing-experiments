# Fampa–Lee AR(1) solid-angle infimum

`check.py` uses SymPy integers and rationals to reproduce the AR(1) precision/covariance identities, the facet point (b^*=(1,0,-\rho^2)^\top), the Van Oosterom–Strackee trihedral-angle calculation, and the Gaussian-orthant denominator formula from Section 8.1 of Fampa and Lee, arXiv:2610.12665v1.  It then evaluates the exact boundary family

\[
\xi_\varepsilon=(1+\rho^2+\varepsilon,\rho+\varepsilon,\varepsilon),
\qquad b_\varepsilon=A_\rho\xi_\varepsilon,
\]

for (ho=1/5,1/2,4/5) and decreasing positive rational (arepsilon).  The cone angle is obtained from an exact rational cross-section hull and exact trihedral triangulation.  Decimal lines are high-precision renderings of exact SymPy expressions and are diagnostic only.
