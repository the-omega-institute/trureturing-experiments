# Gutzeit–Shaban–Yeats–Zalel Covtree witness-bound check

`check.py` uses only exact integer and `Fraction` arithmetic. It reproduces the
height-two family count from the source arXiv:2605.00622v1 and checks
the proposed Section 6 bound on the family

\[
Q_r=L_r\sqcup L_r\sqcup L_{r+1},\qquad n=2r+2,
\]

where `L_s` is `s` incomparable minimal elements below one maximum. The four
size-`n` downset types are enumerated symbolically by component choices. The
forcing certificate uses the `A_n`, `L_r\sqcup L_r`, and
`L_{r+1}\sqcup A_r` downsets: the `A_n` copy supplies `n` global minima, and
overlapping predecessor sets would create a forbidden two-maximum type after
the required unused minima are added to reach size `n`. Thus the three stars
are disjoint and form a downset of size `3r+4`.
