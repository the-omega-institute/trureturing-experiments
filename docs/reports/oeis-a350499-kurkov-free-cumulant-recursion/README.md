# OEIS A350499: Kurkov free-cumulant recurrence

`check.py` performs an exact symbolic check through index 8. It constructs the
compositional inverse of
\[
M(t)=t+\sum_{j\ge1}m_jt^{j+1},
\]
forms the shifted reciprocal `t/M⁻¹(t)`, and compares every coefficient with
Kurkov's recurrence
\[
R(1,k)=m_k,\qquad R(n,k)=R(n-1,k+1)-\sum_{j=1}^{n-1}R(j,k)R(n-j,1).
\]
It also reproduces the five explicit cumulants printed in OEIS A350499. All
arithmetic is symbolic over SymPy expressions; no floating point values are
used. Run with `python3 check.py`; success ends with `ALL_CHECKS_OK`.
