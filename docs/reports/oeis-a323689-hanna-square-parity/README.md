# Hanna's Gaussian Lambert-series parity checks

The exact checks concern Paul D. Hanna's OEIS [A323689](https://oeis.org/A323689)
and [A324300](https://oeis.org/A324300). A323689 prescribes a constant term of 1;
A324300 is locally finite without that replacement.

Run with Python 3 from any directory:

```sh
python3 docs/reports/oeis-a323689-hanna-square-parity/check.py
```

`check.py` uses only Python integers and standard-library functions. Gaussian
integers are pairs of integers. The coefficient formula obtained by multiplying
the numerator's binomial expansion by the negative-binomial denominator expansion
is checked against the symmetric Delannoy formula. Polynomial division independently
checks the kernel for all m,k in 0..24 and the source generating series through degree
80. Delannoy numbers are also checked by their lattice-path recurrence.

`output.txt` records agreement with all 10,001 published A323689 terms (0..10000)
and all 10,202 published A324300 terms (0..10201), with no coefficient, integrality
or parity mismatch. The respective parity targets are N square and N+1 square.
The final line is `ALL_CHECKS_OK`. These finite checks do not establish either
unbounded conjecture.

The input tables `b323689.txt` and `b324300.txt` are Paul D. Hanna's tables from
<https://oeis.org/A323689/b323689.txt> and
<https://oeis.org/A324300/b324300.txt>, respectively, retrieved at
2026-10-25T02:47:07Z. Their SHA256 values appear in `output.txt` and identify the
input bytes used by the program. OEIS material is subject to the
[OEIS End-User License Agreement](https://oeis.org/LICENSE); no new license grant
for these tables is asserted here.
