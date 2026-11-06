# Regionally proximal backward inclusion

Source: Angelina Blahodatna, Lauren Detmold, Daniel Glasscock and Anh N. Le,
*The regionally proximal relation for commutative semigroup actions*,
[arXiv:2610.10397v1](https://arxiv.org/abs/2610.10397v1), Question 7.1.

For a continuous semigroup action on a compact Hausdorff space, RP is the
intersection over open diagonal neighborhoods of closures of unions of diagonal
preimages. RPMinus uses diagonal images instead. The proposed negative answer to
RPMinus being contained in RP uses positive iterates of

\[
f_t(x)=x+t x(1-x^2),\qquad -1\le x\le1,\quad 0<t<1/2.
\]

The proposed relations are RP = R and RPMinus = I², where
R = [-1,0]² ∪ [0,1]². The mechanism is attraction to opposite endpoints in
forward time and attraction of every interior point to zero in backward time.
Time reversal exchanges the relations. The product of this map with its inverse
has incomparable relations on the connected square.

Run with Python 3 and SymPy:

```sh
python3 docs/reports/regionally-proximal-backward-inclusion/check.py
```

The stored standard output is `check.out`. The program uses exact rational and
symbolic arithmetic throughout, and checks:

- Source fidelity: Example 3.3 on the discrete set {0,1}, with the constant map
  to 1, gives RP = X² and RPMinus = {(1,1)}. Including an identity instead gives
  diagonal RPMinus, so the positive-iterate convention matters in this example.
- Eight polynomial identities certify the algebra behind strict monotonicity,
  the three fixed points, half-interval drift, oddness, and geometric bounds.
  Their factor signs apply to every real parameter 0 < t < 1/2.
- Twelve rational orbits for five steps, nine inverse root brackets of width at
  most 2^(-40), and twelve geometric witness bounds, with largest exponent 512.
- The two explicit incomparability witnesses in the proposed product relation
  table: ((-1,0),(1,0)) and ((0,-1),(0,1)).

The inverse bound follows from the checked identity

\[
f_t(z)-(1+t(1-r^2))z=t z(r-z)(r+z)\ge0
\quad(0\le z\le r<1).
\]

Consequently inverse iterates of y with |y| ≤ r contract by at most
1/(1+t(1-r²)) at each step. The forward bound is

\[
(1-ta)(1-x)-(1-f_t(x))
=t(1-x)((x-a)+x^2)\ge0
\quad(0<a\le x<1).
\]

The program validates these universal algebraic identities and finite rational
inequalities. It does not certify the topological definitions, limiting
arguments, inverse existence, or full RP formulas. Those require the mathematical
argument and subsequent Lean verification. No finite grid is claimed to prove a
statement about all real parameters or all open neighborhoods.

The check program is original; no third-party proof or program source was
copied. SymPy supplies exact polynomial arithmetic.
