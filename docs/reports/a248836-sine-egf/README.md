# OEIS A248836: sine EGF check

`check.py` uses only integer arithmetic and `fractions.Fraction`. It exhaustively
enumerates the bounded arrays for `n = 0, ..., 9`, reproducing the first ten
values in the OEIS b-file. It then constructs the formal coefficients of
`(1 - sin(x))⁻²` through exact series inversion and multiplication and compares
them with `a(n)/n!` for the same range. The output is stored in `check.out`.
