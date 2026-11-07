# OEIS A347205 finite popcount rows

`check.py` evaluates the scalar recurrence from the OEIS A347205 NAME and the proposed finite row recursion using exact Python integers. It reproduces the first sixteen published terms, checks row length `popcount(n)+1`, checks `a(n)` against the last row entry, and checks `a(2n)` against the row sum for every `0 <= n <= 16384`. It also rechecks the scalar even recurrence on all tested positive odd half-indices. `check.out` is the captured output; the final line is `ALL_CHECKS_OK`.
