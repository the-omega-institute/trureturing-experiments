# Manjunath polyhedral discriminant check

`check.py` is a Python 3 exact-arithmetic reproduction of the Section 7
quartic family in Madhusudan Manjunath, *Equations over Polyhedral Semirings*,
arXiv:2410.15426v1. It expands

\[
  (Y+[2,3])\oplus(Y+[2,4])\oplus(Y+[-1,1])\oplus(Y+[0,1])
\]

using rational interval endpoints, evaluates all 16 monomials in the ordinary
quartic discriminant support, and checks the root interval `[-2,38]`. It also
checks the affine family `P'_j=t+lambda P_j` and the exact endpoint-slope
obstruction: the lower endpoint forces `u=2`, while the upper endpoint forces
`v=1`, so a full-support multiplicity-two interval cannot exist. No floating
point arithmetic is used.

Run with `python3 check.py`; the final line is `ALL_CHECKS_OK`.
