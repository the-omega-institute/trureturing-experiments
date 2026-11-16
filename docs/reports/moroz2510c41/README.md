# Moroz continuity conjecture check

`check.py` is an independent exact-arithmetic check for the constant-
`P = 1` Lüroth specialization of Moroz's Section 4.1 conjecture. It verifies
the explicit jump value at `s_1 = 6/11`, the stated right limit, the nonzero
jump formula for `s_a`, and the positive Lüroth orbit of `r_a` for `1 <= a <=
100`. The implementation uses only Python integers and `fractions.Fraction`;
it uses no floating point arithmetic.
