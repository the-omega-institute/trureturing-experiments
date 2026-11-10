# Shapiro Conjecture 1.2(a): midpoint nonconvergence candidate

Source: Boris Shapiro, *The enigma of locally Toeplitz matrices. I. Discrete Schrödinger operators with complex potential*, arXiv:2611.09301v1, equation (1.1) and Conjecture 1.2(a). The literal indexing is `0 <= i,j < N`, with entry `a[j-i]((i+j)/(2*N))`.

Run with Python 3.12 or newer and SymPy 1.14.0:

```sh
python3 docs/reports/shapiro-2611.09301-conjecture-1.2a-midpoint-nonconvergence/check.py
```

The coefficient functions are `a[0](x)=1-2*x`, `a[1](x)=x-1/4`, `a[-1](x)=x-3/4`, and zero otherwise. Every mathematical check uses integers, rational numbers or exact symbolic expressions. `output.txt` is the program's complete standard output; its final line is `ALL_CHECKS_OK`.

The fidelity check reproduces Section 1.3's eigenvalue formula for a frozen symmetric tridiagonal matrix of size three with symbolic diagonal and adjacent entries. For the proposed counterexample, the program cross-checks dense matrices against an independent characteristic-polynomial recurrence for sizes 1 through 12; checks the shifted rising-factorial determinant formula for all leading sizes in matrices of sizes 1 through 50; and verifies the determinant, nilpotent block, translation matrix and subsequence separation for every `m=1,...,12`. Factorial estimates are checked for `r=1,...,128`. Three nonzero complex affine transformations are checked for sizes 1 through 8.

With `epsilon=1/559872`, exact polynomial coefficient dominance certifies that the matrices of sizes `4*m` in this finite range have no eigenvalues in the closed epsilon disk. Those of sizes `4*m+2` have at least half their eigenvalues exactly zero, counting algebraic multiplicity. Thus the fixed bounded continuous test `max(0,1-|z|/epsilon)` has integral zero on the checked `4*m` subsequence and at least one half on the checked `4*m+2` subsequence. The weaker determinant-product estimate needed by the infinite-family outline is also checked. No floating point eigenvalue calculation is used.

These are finite exact checks supporting a preregistration. They do not prove the assertion for all sizes, establish literature novelty, or supply a Lean proof. The proposed general argument uses block determinants, a polynomial differential operator, the rising-factorial continuant identity and an elementary factorial lower bound. No third-party program source is incorporated.
