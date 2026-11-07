# OEIS A359506 endpoint bound

`check.py` implements the OEIS incremental reachable-XOR computation with Python integers. It checks the published initial values, independently computes every exact minimum for `n = 0..512`, verifies the proposed closed formula and the bound `3 E(n) <= 5 (n+1)`, and checks the sharp equality family `n = 3*2^s - 1` for `s = 0..12`. The output is recorded in `output.txt`; no floating-point arithmetic is used.
