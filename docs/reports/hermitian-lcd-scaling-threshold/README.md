# Hermitian LCD scaling thresholds

The target is the exact-threshold question in Section 8 of Rahim
Rahmati-Asghar, *Geometric Obstructions and Sharp Genus Dependence in LCD
Scalings*, arXiv:2611.04863v1. For the dimension-three code on
`y^r+y=x^(r+1)` over `GF(r^2)`, the proposed threshold is
`r^3+r^2-2r-2` for every power of two `r >= 4`.

Run `python3 check.py`. The program uses integer bit-polynomial arithmetic
in finite fields and verifies its irreducible moduli by polynomial trial
division. There are no external dependencies or floating-point operations.
Its output is in `check.out`.

The fidelity check reproduces Theorem 5.2: on the elliptic curve over
`GF(4)`, the multiplier `1+x^3+xy` has pole order six and gives the exact
Gram matrix printed in the paper, with determinant one. All lower budgets
are excluded by vanishing moments on a complete monomial basis.

For `r=4` and `r=8`, the program checks every rational point, the indicator
polynomials at three selected points, both moment cancellations, the full
expanded multiplier, its exact weighted pole order, nonzero evaluations,
the information-set determinant, and the scaled Gram determinant.
It obtains thresholds 70 and 558, respectively.

For the lower bound, the three linear functionals `sum_D(f*u)` for
`u=1,x,y` vanish on every basis vector of `L((m0-1)Q)`. By linearity they
vanish on the entire space, including every smaller pole budget. In
characteristic two the Gram diagonals are their squares. A separate
symbolic determinant expansion checks that a symmetric three-by-three
zero-diagonal matrix has determinant zero. This covers every multiplier
in these spaces, rather than sampling coefficient tuples.

The program also verifies the source's residue-sum identity on the entire
basis of `L((n+2g-2)Q)` for both fields. The tests are exact checks of these
two cases. The proposed theorem for all powers of two requires a general
proof; no Lean proof or claim of kernel verification is supplied here.

The check program is original code; its mathematical inputs are the source
paper and the explicit three-point construction described above.
