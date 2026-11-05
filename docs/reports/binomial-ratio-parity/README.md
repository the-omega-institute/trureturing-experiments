# Binomial-ratio parity

`check.py` uses Python 3 standard-library integers and `fractions.Fraction`.
Run `python3 check.py`; `output.txt` contains its output. A successful run
ends with `ALL_CHECKS_OK` and exits 0.

The source is [OEIS A295462](https://oeis.org/A295462), internal record #19,
with Peter Bala's parity conjecture and linked divisibility note. The twelve
published terms provide a fidelity check; the two divisibilities are also
checked through index 4095. OEIS data are subject to the
[OEIS End-User License Agreement](https://oeis.org/LICENSE).

The program checks the conjectured parity characterization against exact
factorial ratios through index 4095, with separate factorial/Fraction
calculation through index 256. It checks factorial valuations, digit sums,
spacing, and first-bad-pair witnesses through index 65535. The 30-cell
Landau table and three rational interval bounds use exact arithmetic.
These are finite checks and do not certify the unbounded conjecture.
