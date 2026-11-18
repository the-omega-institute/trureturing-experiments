# Kushkuley’s pointwise character inequalities

`python3 docs/reports/kushkule2510/check.py` runs with Python 3.12 or later
and uses the standard library. Its output is in `check.out`.

The quaternion multiplication table reproduces the $5/8$ commuting
probability in Example 1 of *Some Remarks on Commuting Probability*.
For odd $m$ from 5 through 21 the program counts every ordered pair in
the dihedral group of order $2m$. Rational polynomial arithmetic modulo
the $m$th cyclotomic polynomial verifies the character calculations,
including the Frobenius count and equation (2.4) for $1+\chi$ and its
extensions by one and three copies of the regular character. Rational
bounds certify strictness of (2.5) at every element and a negative
value at a specified rotation. Equation (2.6) is (2.5) at the identity.
The output gives lower bounds for the normalized gaps of $1+\chi$ and
upper bounds for its witness values. The exact calculation of the
unnormalized gap verifies its invariance under the regular extensions.

The checked sample covers nine group orders. The preregistration's
symbolic argument covers every odd $m\geq 5$ and every positive integer
number of copies of the regular character.

[arXiv:2510.04458v1](https://arxiv.org/html/2510.04458v1)
