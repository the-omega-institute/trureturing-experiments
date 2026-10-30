# OEIS A401145 local-maximum check

`check.py` enumerates unordered partitions into distinct positive parts by exact integer recursion, tests `n + product(parts)` with SymPy's exact primality predicate, and reproduces the OEIS values at `n = 7, 8, 9`. It checks the proposed counterexample at `n = 8`: `a(8) = 1` while `a(7) = a(9) = 3`, and the counted partition `(5, 3)` gives the prime `23`.

The saved output is `check.out`.
