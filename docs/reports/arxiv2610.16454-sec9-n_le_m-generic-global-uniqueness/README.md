# Refutation witness for the N <= m extension of Proposition 4.16

This report checks, with exact rational and symbolic arithmetic, the two-site/two-mode family proposed against the open case stated in Section 9 of Juan Antonio Vega Coso, *Axioms and internal rigidity of spectral duality structures*, arXiv:2610.16454v1. It uses

- `A = [[1,2],[2,1]]`, `pi = (1/2,1/2)`, `beta = 1`, and spectral functions `f(gamma) = (1,gamma)`;
- `D = [[a,-a],[b,-b]]`, `B = A + D`, with `0<a<1`, `0<b<1`, `a/2<b<2a`, and `a != b`;
- `q = ((2b-a)/(a+b),(2a-b)/(a+b))`.

The script proves symbolically that `D*pi=0`, `D*q=(b-a)A*q`, `C(pi,gamma)=1/2`, and `C(q,gamma)=1/(2+b-a)`. The exact point `a=1/2`, `b=3/4` has `q=(4/5,1/5)`, `C(pi)=1/2`, and `C(q)=4/9`. It also checks positivity, H1, A-column independence, and failure of B2 by exact determinants. The strict inequalities define a nonempty open subset of the full linear space `D*pi=0`, so the N>m generic uniqueness conclusion cannot extend to N<=m.
