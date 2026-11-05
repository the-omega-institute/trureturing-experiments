# Minimal-upper-bound reachability checks

Source: Adam Humeniuk, Gabriel Jarry-Bolduc, Patrick Pascua and Nejaunie
Williams, *An algorithm to exactly compute minimal upper bounds in the Loewner
order*, arXiv:2606.18173v1, Algorithm 3.4, Example 3.9 and Remark 3.10.

Run with Python 3 and SymPy:

```sh
python3 check.py
```

`check.out.txt` contains the output. All calculations use exact integers or
Gaussian rationals. Positive semidefiniteness is checked by every principal
minor. Each rank-one coefficient has both a feasibility check and the exact
quadratic-form identity `u*(P-s ee*)u = q-s q^2`, certifying maximality.

The program reproduces Example 3.9's two updates, coefficients and terminal
matrix. It tests the proposed direction `e=(B-M)v`, with `v` in a target contact
kernel, on 1,076 prescribed-target instances and verifies feasibility, agreement
with Algorithm 3.4's coefficient, terminality, and the residual rank drop.

The finite grid consists of all scalar families of sizes 1, 2 and 3 with entries
in -2 through 2 and four starting offsets, together with the declared real and
complex contact/residual grids in dimensions 2 and 3. In the matrix grids,
`H_i=T* diag(h_i) T` has its i-th diagonal weight zero; all other weights are
positive. The targets therefore have spanning contact kernels. Residuals are
`R diag(d) R*`, including singular residuals and zero-length runs. Every declared
grid case is checked; this is not an enumeration of all Hermitian matrices.

An ancillary exact instance shows that the contact-space dimension may increase
from zero to two in one update, contrary to Remark 3.8's exact-one assertion.
This does not refute Remark 3.10's reachability assertion. The tested
endpoint-dependent count is `rank(B0-M)`.

These checks supply finite evidence and symbolic form identities. They do not
prove the universal reachability claim; its proposed settlement uses a general
finite-dimensional argument rather than finite positive enumeration.
