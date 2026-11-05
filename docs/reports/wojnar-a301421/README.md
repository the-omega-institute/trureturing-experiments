# Wojnar A301421: exact finite checks

Run with Python 3, using the standard library only:

```sh
python3 docs/reports/wojnar-a301421/check.py
```

The program reproduces the first three six-data normalized Girard–Waring polynomials and the first seven positive-coefficient sums in Table GW.n=6 of Wojnar, Wojnar and Brin, arXiv:1706.08381v2, equation (2) and p. 24. It also reproduces all 22 DATA terms in OEIS A301421 revision 31. The sequence has offset one.

For every degree 1 through 24, two independent exact constructions produce the same monomial coefficient dictionary: partition enumeration using the source factorial formula, and Newton recursion using signed binomial-normalized elementary coefficients. Positive coefficients are selected after monomials have been aggregated. Their masses agree with the coefficients of the published rational generating function. Exact evaluation on four six-data sets checks the source power-mean interpretation, including repeated, zero, negative and rational data. Arithmetic uses integers and fractions.Fraction throughout.

The stored output ends with ALL_CHECKS_OK. These finite checks support source fidelity and test the expected positive settlement; they establish no all-degree identity or Lean theorem. The target is only the six-data GF assertion. The separate ratio assertion is outside these checks.

Sources:

- https://arxiv.org/abs/1706.08381v2
- https://oeis.org/A301421/internal
- https://oeis.org/A301417/internal

Original check code is offered under the MIT license; source mathematics is cited, and the report redistributes no source manuscript.
