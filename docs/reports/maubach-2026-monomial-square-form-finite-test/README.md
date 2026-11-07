# Maubach Problem 9 exact check

This report checks the negative finite-test construction for Stefan Maubach,
*Universality of monomial square forms*, arXiv:2611.03057v1, Problem 9.

For `F_k(x) = sum_{j=1}^k (x_1 ... x_j)^2`, the monomials are distinct,
nonconstant, squarefree, and have coefficient one.  Setting the first `r`
variables to `1` and the rest to `0` gives `F_k = r` for every `1 <= r <= k`.

If `F_k(x) <= 2^(2t)` and the first `j` prefix product is positive, each
factor has a dyadic bin `2^b <= |x| < 2^(b+1)` with sum of bin indices at most
`t`.  There are at most `(t+1)^j` bin vectors, and each bin contains at most
`2^t` integer tuples.  Summing over `j <= k` gives the upper bound
`k * 2^t * (t+1)^k` on the number of distinct positive values at most
`2^(2t)`.  With `t = 4 k^2`, the exact inequalities
`k <= 2^(k-1)` and `4 k^2 + 1 <= 2^(3k)` give a bound at most
`2^(8 k^2 - 1)`, less than half of the available positive integers.  For
`k >= 4` this also leaves a missed value greater than `k`.

`check.py` uses only exact integer arithmetic (and `fractions.Fraction`),
enumerates negative, zero, and positive inputs for `k <= 4`, checks the
specializations and the exact power inequalities, and prints `ALL_CHECKS_OK`.
Run it from this directory with:

```text
python3 check.py | tee output.txt
```
