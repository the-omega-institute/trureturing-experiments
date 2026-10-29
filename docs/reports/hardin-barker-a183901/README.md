# Hardin–Barker OEIS A183901: exact prospective checks

Source: R. H. Hardin, [OEIS A183901](https://oeis.org/A183901/internal),
revision 8 (5 April 2018), including Colin Barker's two FORMULA conjectures.
The sequence definition and data are attributed to Hardin; the generating
function and recurrence are attributed to Barker. The 37 reference terms come
from [Hardin's b-file](https://oeis.org/A183901/b183901.txt). OEIS content is
available under [CC BY-SA 4.0](https://oeis.org/LICENSE).

Run with Python 3, using only its standard library:

```sh
python3 docs/reports/hardin-barker-a183901/check.py
```

The stored result is `output.txt`. The check uses integers throughout and exits
0 only after printing `ALL_CHECKS_OK`.

The program checks all 78,125 multiplicity vectors capped at four, computes the
exact finite numerator identity, and compares their cap-fiber counts with all
37 source terms. An independent dynamic program selects three distinct
individual occurrences after excluding the target occurrence. It exhausts all
sorted arrangements for n = 1 through 7 and reproduces the source counts
1, 1, 136, 664, 1720, 3491, 6263. It also checks the claimed recurrence through
n = 300 and its nonzero residual −33 at the excluded index n = 10.

For cap vector d, write s = sum(d) and r = number of entries equal to four.
The polynomial check is

    sum_valid_d x^s (1-x)^(7-r) = (1-x)^7 + x^3 Q(x),

where Q is Barker's complete numerator including the initial x. The program
checks coefficients through degree 28, including all zero high coefficients.

These are exact experimental checks, not a Lean proof. A proof for every
length requires the sorted-arrangement/multiplicity bijection, cap invariance,
the independent residual-multiplicity fiber bijection and formal-series
algebra. The finite identity is a helper to that general argument. Bounded
recurrence checks alone give no statement about untested indices.
