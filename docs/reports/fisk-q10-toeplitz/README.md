# Fisk Question 10: exact Toeplitz minor checks

Sources: Steve Fisk, [Questions about determinants and polynomials](https://arxiv.org/abs/0808.1850v1), Question 10, and [Polynomials, roots, and interlacing](https://arxiv.org/abs/math/0612833v2), Section 21.10, Conjecture 21.84. The former uses upper-half-plane nonvanishing for the input and closed-left-half-plane roots for a stable minor.

Run with Python 3 and SymPy 1.14.0:

```sh
python3 docs/reports/fisk-q10-toeplitz/check.py
```

`output.txt` contains the exact run output. The final line is `ALL_CHECKS_OK` and the exit code is 0. All evidence uses integers, rational numbers and exact algebraic expressions; no numerical root finding occurs.

Fidelity: all 330 minors of the first four rows and seven columns of the book's displayed Toeplitz matrix for `(1+x+y)^4` equal nonnegative integer multiples of nonnegative integer powers of `1+x`. There are 264 nonzero minors. This checks the displayed example on that finite block; the book's assertion about every finite minor is broader.

Candidate: for `t>=4`, set `f_t=(t*x+y+1)*(x+y+t)`. The code checks symbolic identities for the six positive coefficients and the imaginary parts of both affine factors, then constructs the literal coefficient minor with rows `[0,1,2,3,4]` and columns `[1,2,3,4,5]`. Independent permutation and SymPy determinant algorithms agree, as does the continuant recurrence. Its factor `Q_t=(t^2-t+1)*x^2-(t^2-4*t+1)*x+(t^2-t+1)` has the exact root `z_t=(C+I*sqrt(3)*S)/(2*A)`, where `A=t^2-t+1`, `C=t^2-4*t+1`, `S=t^2-1`. The program verifies the root identity symbolically and checks that each of `A,C,S` has positive constant term and nonnegative coefficients after `t=m+4`, so `Re(z_t)>0` for `m>=0`. The same root is rechecked for every integer `4<=t<=40`.

At `t=4`, the determinant is `5*(x+1)*(21*x^2+33*x+21)*(13*x^2-x+13)`, and the exact root is `(1+15*I*sqrt(3))/26`, with real part `1/26`.

These computations support a proposed refutation. They are not a Lean proof, a complete classification of minors, or a literature novelty claim.
