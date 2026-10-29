# OEIS A088500: rising-factorial integral

Source: [OEIS A088500](https://oeis.org/A088500), revision 24 (12 January 2026).
The NAME defines the exponential generating function
$1/(1+2\log(1-z))$. The FORMULA section gives
$a_n=\sum_{k=0}^n {n\brack k}k!2^k$ and Velin Yanev's conjecture
$a_n=\int_0^\infty (2x)^{\overline n}e^{-x}\,dx$.
Here the Stirling numbers are unsigned and the empty product is 1.

Run with Python 3, using only its standard library:

```sh
python3 docs/reports/yanev-a088500-rising-factorial/check.py
```

`check.out` contains the output. All arithmetic uses integers or exact rational
numbers. For n = 0 through 50, three separate constructions agree: inversion of
the defining EGF denominator, the weighted Stirling sum, and the constant term
of the exact polynomial Q satisfying Q-Q' = product(2x+j, j=0..n-1).
The last construction is the boundary value of the exact antiderivative
-exp(-x)Q(x). The program verifies every coefficient of that identity.
It also reproduces all 20 initial values printed in the source and its
Gutkovskiy recurrence, including the n=0 empty-product case.

These are finite checks. The program does not establish the general identity,
integrability, or the polynomial-times-exponential limit at infinity. Those
analytic facts and the universal quantifier remain obligations of a later
proof. No numerical quadrature is used.

Related method: Sela Fried, [Proof of a conjecture stated in A007840](https://oeis.org/A007840/a007840.pdf),
Theorem 1, treats the coefficient-1 version. The check program is original
code; its mathematical input and initial values are attributed above.
