# Krishna's C*-algebraic Sendov conjectures: a branch-exchange counterexample

`check.py` checks the cubic counterexample over A = C([0,6], ℂ) to Conjectures 2.4 and 2.5 of arXiv:2203.06916v1.
The roots are a₁ = 1/2, a₂ = −1/2 and a₃ = c, where c is the closed rational polygon through −9/10, −1/10+4i/5, 1/10+4i/5, 1/10+9i/10, −1/10+9i/10, −1/10+4i/5, −9/10.

The script checks, on a grid of 60001 points:
- the two critical roots b± = (c ± √(c²+3/4))/3, continued along the loop, solve p′(b) = 3b² − 2cb − 1/4 = 0;
- the two critical roots are positive barycentric combinations of the roots and stay inside the unit disc;
- the discriminant and the root gaps stay away from 0;
- the branches exchange between t = 0 and t = 6;
- |b₋(0) − 1/2| = |b₊(6) − 1/2| = (12+√39)/15 > 1.

It also checks two exact facts with rational arithmetic: the vertex moduli are at most √82/10, and the imaginary-axis crossings 4i/5 and 9i/10 straddle the branch point i√3/2.

Run: `python3 check.py` (exit 0, final line `ALL_OK`).
