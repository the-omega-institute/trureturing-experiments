# Mathar's order-four recurrence for OEIS A124431

`check.py` evaluates the literal floor-binomial sum using Python integers. It
compares all 25 displayed OEIS terms and independently counts peakless Motzkin
paths with three horizontal colors at height zero and two at positive height.
The path computation reproduces Deutsch's example a(3) = 27 + 2 = 29.

The program checks both Mathar recurrence transcriptions for n = 4 through 250,
the parity split for n = 0 through 100, and the constant-term identities through
degree 60. It checks the central-binomial convolution through degree 100 and the
proposed differential equation through degree 250. No floating point arithmetic,
external Python packages, network access, or generating-function assumptions are
used by the program. These finite checks do not prove the unbounded recurrence.

Run from this directory:

```sh
python3 check.py
```

The successful output is `check.out`; its final line is `ALL_CHECKS_OK`.

Source: [OEIS A124431](https://oeis.org/A124431), revision 41. The sequence is due
to Paul D. Hanna; the recurrence is conjectured by R. J. Mathar in COMMENT and
FORMULA. The generating functions already recorded by Maksym Voznyy and Paul
Barry are prior work. This report tests a derivation from the literal defining
sum and does not assume either generating function. The check program is original
code under the repository's Apache-2.0 license.
