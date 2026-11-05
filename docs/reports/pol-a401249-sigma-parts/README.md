# A401249: minimum subparts per symmetric-sigma part

Run from the experiments repository root:

```sh
python3 docs/reports/pol-a401249-sigma-parts/check.py
```

Python 3 standard library only; all arithmetic is exact integer arithmetic.
`check.out` is the program output and ends with `ALL_CHECKS_OK`.

The geometry is computed from the triangular-offset columns of
[A237048](https://oeis.org/A237048), their alternating partial sums
[A249223](https://oeis.org/A249223), and the reflected width rows
[A262045](https://oeis.org/A262045). Parts are positive contiguous runs;
subparts are connected components of each width-one layer, following
[A279387](https://oeis.org/A279387). The arithmetic side separately enumerates
positive divisors and splits at consecutive ratios greater than two.

The program reproduces the n=70 example and all 105 published DATA terms of
[A401249 revision #22](https://oeis.org/A401249), plus the layer examples at
n=12, 15, 18, 45, 90, 360 and the per-part examples in
[A280940](https://oeis.org/A280940). For every 1 <= n <= 10000, the ordered
per-part subpart counts equal the ordered per-block odd-divisor counts.
This covers 32598 parts and 50292 subparts. It also checks the published
half-open divisor-window width identity and the local positive-variation
count used in the proposed general proof.

These are finite checks, not a universal proof. The layer calculation does
not define geometric counts using the conjectured divisor-block equality.
The general conjecture and its stronger A280940 formulation remain the
targets of a later proof stage.

Sources: OEIS entries above and Höft's
[width identity](https://oeis.org/A249223/a249223.pdf) and
[part/block correspondence](https://oeis.org/A237270/a237270_4.pdf).
Sequence data and short source examples are subject to the
[OEIS End-User License Agreement](https://oeis.org/LICENSE).
The check program is original code under the repository's Apache-2.0 license.
