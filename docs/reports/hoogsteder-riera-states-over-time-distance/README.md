# Hoogsteder–Riera stability-factor check

This report checks the proposed improvement in question (3) of Section 5 of
Matt Hoogsteder-Riera, *Membership in the cone of quantum states over time is
NP-hard*, arXiv:2610.34972v1 (26 October 2026).

`check.py` uses only integers, `fractions.Fraction`/SymPy rationals, and exact
algebraic eigenvalues. It first reproduces the source pairing formula on an
exact qubit unitary channel. It then enumerates every nonzero vector in
`{-1,0,1}^{d r}` for `(d,r)=(2,2),(3,2),(2,3)`, forms the rank-one Hermitian
matrix `Z = z z†`, and checks the proposed universal estimate

`||{rho tensor I, Z^T_A}/2||_1 <= sqrt(d) ||Z||_1`

by exact squaring. Finally, an exact isometric channel with defect `eta >= d`
checks the proof's separate high-defect branch
`||Q||_1 <= d sqrt(eta)`.

Run:

```text
python3 check.py
```

The final line is `ALL_CHECKS_OK`.
