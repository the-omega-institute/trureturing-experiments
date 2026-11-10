# Generalized-wheel squarefree-power check

This report independently evaluates the maximum cardinality of a q-free face
(the largest vertex subset whose induced matching number is below q) for the
cycle joined to either a clique or an independent hub set.  The exact integer
matching recurrence checks every `3 ≤ n ≤ 8`, `1 ≤ m ≤ 4`, and every valid
power q, including the zero-generator case q > ν.  It also reproduces the
source's q = 1 dimensions.  The script uses only integer bitmasks and counts;
there is no floating-point evidence.

Run `python3 check.py`; the final output line is `ALL_CHECKS_OK`.
