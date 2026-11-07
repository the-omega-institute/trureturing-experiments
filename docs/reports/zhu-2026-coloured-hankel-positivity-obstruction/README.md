# Coloured Hankel positivity obstruction

This report checks the coloured functional-digraph EGF from Bao-Xuan Zhu,
*Coefficientwise Hankel-total positivity from a universal cumulant lifting*,
arXiv:2611.01172v1.  It uses exact rational formal power-series arithmetic.

Run:

```sh
python3 check.py
```

The script verifies the tree-function equation through degree 14 and checks,
for two and three independent colours and every `1 <= n <= 6`, that the
coefficient of `y_1 y_2 q_1^(n+1) q_2^(n+1)` in
`M_n M_(n+2) - M_(n+1)^2` is `-2 (n!)^2`.  This determinant is the Hankel
minor with rows `{0,1}` and columns `{n,n+1}`.  The negative terms come from
the two ordered cross-products in the square of the degree-one component
intensities.  The positive product has degree at most `n` in one factor and
therefore cannot supply both target exponents.

The source's established forms remain distinct: fixing all `q_l >= 0` gives
coefficientwise positivity in the `y_l`; setting all colours equal gives the
one-colour question; and the numerical moment statement is weaker than joint
coefficientwise positivity.
