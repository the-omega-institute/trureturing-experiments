# Bala A337389 all-column diagonal check

Run `python3 check.py`; the expected output is `output.txt`, ending in `ALL_CHECKS_OK`.
The program uses Python's standard library with exact integers and `fractions.Fraction`.

Source: [OEIS A337389](https://oeis.org/A337389), revision #40, May 05 2022.
The array is attributed to Seiichi Manyama; the final FORMULA conjecture is attributed
to Peter Bala, May 03 2022. Source text and displayed numerical data are available
under the [OEIS End-User License Agreement](https://oeis.org/LICENSE).

For n=0..30 and each of 36 distinct rational parameters (k=0..30, 101, 1000000,
-1, -4, and 1/2), the program compares three independently specified expressions:
the source's two normalized formal square roots, the conjectured coefficient
of the n-th power of its rational series, and the source's finite binomial sum.
It reproduces the 36 displayed array entries and all 55 published prefix entries,
and checks 1044 instances of the source recurrence. All 1116 identity pairs agree.

The square-root calculation determines each coefficient recursively from a square;
the diagonal calculation repeatedly multiplies the independently specified rational
series. Each truncated root and inverse is checked by multiplication. The finite
binomial sum supplies a third computation. No floating point is used.

These results provide source fidelity and finite confirmation. They do not prove
the identity for unbounded n or k, settle its literature status, or constitute a
Lean proof. The negative and fractional parameters test the proposed algebraic
mechanism and are outside the OEIS assertion's nonnegative integer parameter range.
