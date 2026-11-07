# OEIS A183683 exact check

`check.py` enumerates all 32 binary rows, tests the four consecutive 2-by-2 integer determinants exactly, and counts labelled height `n+1` arrays by dynamic programming. It reproduces the first ten OEIS terms, verifies the four-class quotient matrix, checks the recurrence through source index 200, and checks the first coefficients of the formal-series identity. No floating-point arithmetic is used.

Run with Python 3.12 or newer:

    python3 check.py

The expected final line is `ALL_CHECKS_OK`.
