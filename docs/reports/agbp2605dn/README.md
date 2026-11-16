# agbp2605dn exact arithmetic check

`check.py` uses only Python integers, `fractions.Fraction`, and SymPy exact
radicals. It verifies the capsule radius equations, the exact perimeter and
shell-volume identities, and the rational separation used by the refutation:
the capsule Rayleigh lower bound is strictly greater than `k/3` for `k ≥ 729`,
whereas the comparison annulus has a trial quotient below `108 < 243`.

Run with `python3 check.py`. The stored output is `check.out`; its final line
is `ALL_CHECKS_OK`.
