# Moroz continuity conjecture check

`check.py` is an independent exact-arithmetic check for the constant-
`P = 1` Lüroth specialization of Moroz's Section 4.1 conjecture. It verifies
the explicit jump value at `s_1 = 6/11`, the stated right limit, the nonzero
jump formula for `s_a`, and the positive Lüroth orbit of `r_a` for `1 <= a <=
100`. The implementation uses only Python integers and `fractions.Fraction`;
it uses no floating point arithmetic.

These finite checks verify the candidate's algebra; they do not verify a limit,
continuity, membership in the exceptional set, or the full source conjecture.
The literal conjecture is already contradicted in the literature at the endpoint
1: Ahn, *Baire-category transfer between Engel and Pierce expansions*,
arXiv:2610.00946v1, Theorem 3.5(i) proves continuity there, while the paragraph
before that theorem identifies the exceptional set as the rational numbers in
(0,1]. No new preregistration of that literal conjunction is proposed.
