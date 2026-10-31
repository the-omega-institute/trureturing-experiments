# OEIS A206982 exact recurrence check

`check.py` implements the literal four-direction avoidance definition for `n × 3`
binary arrays, and an independent transfer calculation on the 36 legal row pairs.
It checks Hardin's published values through `n = 7`, the order-five recurrence for
`8 ≤ n ≤ 40`, and the coefficient identity for Barker's rational generating function
through degree 40.  The arithmetic is integer-only.

Run with Python 3:

```text
python3 check.py
```

The final line is `ALL_CHECKS_OK`.
