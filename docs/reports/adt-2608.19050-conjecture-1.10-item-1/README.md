# Candidate check: arXiv:2608.19050, Conjecture 1.10(1)

`check.py` is a Python 3 exact-arithmetic check. SymPy verifies the logarithmic Bellman integral identity and the rescaled variational identity symbolically; integer and `Fraction` checks verify the alternating permutations and the two-block decomposition through `n = 100`. No floating-point value is used as evidence.

Run with:

```text
python3 check.py
```

The stored output ends in `ALL_CHECKS_OK`.

As a separate diagnostic reading, a double precision Newton iteration of the
published recurrence through `n = 100000` gives
`b_n / n^2 = 0.4883419588484151`; direct minimization of the variational
function gives `0.4883197421891127`. These decimals are recorded only as a
sanity reading; the check above uses no floating point as evidence.
