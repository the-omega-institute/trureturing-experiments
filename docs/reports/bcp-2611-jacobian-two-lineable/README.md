# Real even-dimensional Jacobian lineability

Source: L. Bernal-González, M. C. Calderón-Moreno and J. A. Prado-Bassas,
*Counterexamples to the Jacobian Conjecture from the viewpoint of lineability*,
arXiv:2611.03572v1, Question 2. The source attributes the displayed map to L. Alpöge.

Run with Python 3 and SymPy:

```sh
python3 docs/reports/bcp-2611-jacobian-two-lineable/check.py
```

`check.out` is the successful output. All polynomial identities and point evaluations
use exact integer or rational coefficients; no floating point is used as evidence.

The program independently reproduces the source's determinant `det DF0 = -2`
and all three printed preimages of `(-1/4, 0, 0)`. It also verifies that the distinct
points `(2, -1/2, 5/4)` and `(-2, 1, 2)` both map to `(0, 1, 0)`.

For dimensions 2, 4, 6, 8, 10 and 12, it checks the symbolic identities
`J^2 = -I`, `det(aI+bJ) = (a^2+b^2)^(N/2)` and
`(aI-bJ)(aI+bJ) = (a^2+b^2)I`, with the standard block-diagonal complex structure.
The determinant is positive at each of 48 nonzero rational coefficient pairs in
each dimension, giving 288 exact comparisons.

For dimensions 4, 6, 8, 10 and 12, it pads the source map by identity coordinates,
checks its determinant and collision, differentiates `(aI+bJ)F` symbolically and
verifies `D((aI+bJ)F) = (aI+bJ)DF`. The two exact determinant identities give
`det D((aI+bJ)F) = -2(a^2+b^2)^(N/2)` by determinant multiplicativity. In dimension
4 the transformed determinant is also computed directly, independently of that
factorization. An invertible derivative at the origin and the pencil entries
recover the two coefficients, supplying a check of linear independence.

These are reproducible algebraic checks for the preregistration. The unbounded
statement for every even dimension and the pointwise construction for arbitrary
members of the source family are a proposed mathematical argument, not a Lean
proof or a consequence of finite sampling alone. The program contains original
code written from the paper's mathematical formulas.
