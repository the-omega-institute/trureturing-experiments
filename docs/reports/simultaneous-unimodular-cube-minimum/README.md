# Simultaneous unimodularity and a fixed-order entry bound

Run with Python 3, from this repository:

```sh
python3 docs/reports/simultaneous-unimodular-cube-minimum/check.py
```

Only the Python standard library is used. All arithmetic is exact: integers,
checked Bareiss divisions, and an independent determinant implementation over
`fractions.Fraction`. The stored output ends with `ALL_CHECKS_OK`.

The mathematical source is Xinyi Li, Shuochen Wang, Zimeng Ye and Alexandru
Zaharescu, *Simultaneous Unimodularity of Positive Integer Matrices and Their
Entrywise Cubes*, [arXiv:2610.14546v1](https://arxiv.org/abs/2610.14546v1),
Section 9, second question. The source cube is entrywise. The program was
written from the paper's displayed mathematics; it does not use its ancillary
program or any third-party implementation.

The fidelity checks reproduce the twenty signed cofactors of the order-ten
seed in equation (2.1), the primitivity and gcd certificates (3.3)-(3.5), and
the reduced Jacobian minor (3.8). Each seed determinant is checked by two
independent exact algorithms. The general bordering identity of Lemma 7.1 is
checked on all 196,830 cases in the program's explicitly bounded grid.

The checks for the proposed route establish all 180 duplicate-row evaluations
used to show that the lower-row cofactor polynomials do not vanish identically
on the homogeneous cone. On all 91,125 small integer matrices specified in the
program, they check Laplace expansion and the degree-one/degree-three cofactor
scaling. On the 6,936 matrices in that grid whose two determinants equal one,
they check the singularity criterion `p_ij = a_ij^2 q_ij` for every row.
The grid includes changed-row replacements over 552 lower blocks. One explicit
replacement and one unit-complement border both preserve the two determinants.

The small algebraic models have entries equal to zero, one, or negative
integers. They verify identities and are not positive order-ten witnesses.
These finite computations do not prove the Section 9 assertion, certify the
circle-method estimates, establish density of the real smooth patch, or
establish unboundedness. The proposed positive settlement uses a general
argument with those separate mathematical ingredients, not positive exhaustive
enumeration. All paper results used as ingredients still require proofs in an
unconditional Lean settlement.
