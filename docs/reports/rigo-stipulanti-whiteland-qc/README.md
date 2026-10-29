# Rigo–Stipulanti–Whiteland Question C check

`check.py` uses integer dynamic programming for scattered-subword counts and
exact `fractions.Fraction` arithmetic. It reproduces the finite inequalities
from Theorem 7.4 for the fixed point of `0 → 01`, `1 → 12`, `2 → 2` at lengths
6 through 15, verifies finite lacunary reconstruction and unbounded abelian
lower bounds for Fibonacci-spaced markers, and checks the proposed product
overlay on finite prefixes for levels 2 through 4. It uses no floating-point
arithmetic.

Run from this directory with `python3 check.py`. The stored output is
`check.out`; its final line is `ALL_CHECKS_OK`. These are finite fidelity and
sanity checks for the preregistration. The universal statement is intended to
be settled by the general lacunary reconstruction, overlay, and synchronized
morphic arguments in the draft, not by exhaustive positive enumeration.
