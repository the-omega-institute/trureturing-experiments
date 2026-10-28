# Han–Pedon higher-shift q-metallic Hankel checks

Source: Guo-Niu Han and Emmanuel Pedon, *Hankel determinants for the
(n+2)-shifted q-metallic numbers*, arXiv:2610.36132v1, Conjecture C.

Run with Python 3; no third-party packages or input files are required:

```sh
python3 check.py
```

`check.out` is the program's complete output. The program uses integers,
including fraction-free Bareiss determinants, and performs these checks:

- Reconstructs the moments from the source quadratic equation and reproduces
  the displayed n=2 period and n=3 antiperiod after Theorem B, together with
  the following block: 72 determinants.
- Checks both transfer matrices, their constant terms, the ordinary transfer
  determinant, the binomial expansion for low jets, and the two Hensel
  branches for high jets, for n=1,...,8 and block counts 1,2,3,7:
  96 finite jet identities. The tested precision is min(2*rho+3,32), and
  omitting the small branch is checked only at degrees strictly below 2*N.
- Checks normal orthogonality and the full Bezout inverse for (n,N) equal to
  (1,4), (1,8), (2,12), (2,24). The determinant at (1,4) is -1;
  at the other three it is +1. Jacobi's formula includes that determinant
  factor. All 252 admissible complementary-minor identities for shift
  offsets 2,...,7 are checked, including empty and zero target minors.
- Checks every residue modulo rho for n=1,2, shifts n+3,n+4, and m=0,...,5:
  240 determinants. Affine Newton polynomials at shift n+3 use the first
  two values and four holdouts. The shift-n+4 check allows degree at most
  three, uses the first four values and two holdouts; the observed degree
  is at most two. These are checks of a finite grid only.

The cycle data are explicit transcriptions of the lower-shift H-fraction
package from arXiv:2502.05993v2; the golden three-step cycle is separate.
No files from the Lean repository are imported or read by the program.
The last output line is `ALL_CHECKS_OK` on success. These computations
corroborate a proposed route; they do not prove universal quasipolynomiality,
minimality of the period, or a universal degree bound.
