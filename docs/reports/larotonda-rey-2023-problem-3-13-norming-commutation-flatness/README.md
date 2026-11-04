# Larotonda–Rey Problem 3.13: exact counterexample check

`check.py` uses only SymPy exact matrices, rational numbers, and radicals. It
reproduces the cyclic trace identity and the Frobenius equality/non-equality
special case used in the source, then checks the proposed two-dimensional
counterexample.

Write a Hermitian matrix as `H=t I+H_0` with

`r(H_0)^2 = Tr(H_0^2)/n`.

For `0<a<1`, the closed, balanced, convex and conjugation-invariant body

`B_a = { iH : t^2 + max(r(H_0)-a,0)^2 <= 1 }`

is the Minkowski sum of the Euclidean unit ball in the central coordinate and
the radius-`a` Euclidean ball in the traceless coordinate. Its gauge is a
nondegenerate, subadditive, positively homogeneous and Ad-invariant Finsler
norm. The script evaluates this gauge on

`K=diag(1,-1)`, `V=i(I+aK)`, `a=1/2`, and
`X=[[0,1],[-1,0]]`.

The unique norming matrix for `V` is `N=iI/2`, so `[N,X_C]=0`. However,
`r((V+s[X,V])_0)^2 = a^2(1+4s^2)` (the exact `q` for this `X,K` is `4`), and
therefore the gauge is strictly larger than its value at `V` for every
nonzero `s`. The script also computes the
positive scale `(raw gauge of iP)^(-1) = 3/(sqrt(7)-1)` needed for the source's
normalization `||iP||=1`; scaling preserves the strict inequality.

Run:

```sh
python3 check.py > check.out
```

The final output line is `ALL_CHECKS_OK`.
