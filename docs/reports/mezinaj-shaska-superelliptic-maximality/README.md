# Mezinaj–Shaska superelliptic maximality

Exact arithmetic checks for the final maximality question in Problem 7.11 of
Jurgen Mezinaj and Tanush Shaska, *Hulls, linear equivalence, and weighted
superelliptic codes*, [arXiv:2608.02850v2](https://arxiv.org/abs/2608.02850v2).

Run with Python 3.8 or later, using only the standard library:

```sh
python3 docs/reports/mezinaj-shaska-superelliptic-maximality/check.py
```

The stored output is `check.out`. Every check uses integers or `fractions.Fraction`.
Fields are coefficient vectors modulo a polynomial whose irreducibility is checked
by exact trial division by all possible monic factors up to half its degree.
Frobenius and multiplicative inverses are checked on every field element.

For r = 3, 5, 7, 9, 11, 13, the program reproduces Theorem 7.10: for every
n >= 2 dividing r+1, the curve y^n = -(x^(r+1)+1) has r+1 rational branch
points, n rational infinity points, genus (n-1)(r-1)/2, and exactly
r^2+1+2rg rational points. The source's evaluation count N excludes branch points.

For the same r, it checks the candidate family n=2 and
f_r(X) = sum(X^(j(r-1)), j=0,...,r). The product identity
(X^r-X)f_r = X^(r^2)-X and the product over the subfield complement are both
checked as polynomial identities. It verifies every value, every root derivative,
the full affine point set, the reciprocal chart, and its overlap with the affine
chart. All r^2-r roots are simple. There are r^2+r affine points and two smooth
rational infinity points. The total counts are:

| r | q | degree | genus | count | Hasse–Weil maximum |
|---|---|---|---|---|---|
| 3 | 9 | 6 | 2 | 14 | 22 |
| 5 | 25 | 20 | 9 | 32 | 116 |
| 7 | 49 | 42 | 20 | 58 | 330 |
| 9 | 81 | 72 | 35 | 92 | 712 |
| 11 | 121 | 110 | 54 | 134 | 1310 |
| 13 | 169 | 156 | 77 | 184 | 2172 |

The base witness over F_9 uses i^2=2 and f=X^6+X^4+X^2+1. The output lists
all evaluations, the six nonzero root derivatives, twelve affine points and two
infinity points. For each tested r, the integer excess over q+1 is r+1, strictly
between 0 and 2r; it cannot equal 2rg for any nonnegative integer genus.

These finite checks support a proposed general refutation. They do not prove
the family for every odd prime power, identify the glued curve with a smooth
projective model in Lean, or settle the classification of all delta=0 pairs.
The polynomial construction is credited to the complementary-factor mechanism
in Borges, [arXiv:1311.4438](https://arxiv.org/abs/1311.4438), Lemma 4.2.
