# Arzhantsev–Zaidenberg Question 2 check

`check.py` uses exact integer arithmetic, `fractions.Fraction` (for the
declared arithmetic layer), and SymPy's exact polynomial ring. It reproduces
the cyclic-weight calculation behind Proposition 5.2/Lemma 5.10, verifies the
classification `ee' = kd + 1`, and checks the proposed mixed-jet obstruction
for the first five `k >= 2` pairs with `d <= 100`, including
`(d,e,e') = (4,3,3)`. The final line of the saved output is
`ALL_CHECKS_OK`. These are independent checks of the arithmetic and finite jet
calculation; the ind-group tangent identification remains a formalization
obligation.
