# Hanna's OEIS A078122 coefficient formula

Source: Paul D. Hanna, [OEIS A078122](https://oeis.org/A078122), second FORMULA, revision 10. The source's displayed matrix and programs use first column one; the NAME's first-row boundary conflicts with lower triangularity.

Run `python3 docs/reports/hanna-a078122/check.py` from the repository root. The program uses Python's arbitrary-precision integers and no external packages. Its output is stored in `check.out`.

The source matrix is built solely from its shift-power recursion. A separate unbounded coin-change dynamic program computes the partition coefficients. The checks reproduce all 45 source data entries through row 8 and the displayed four-by-four cube, compare 144 entries and 32 row sums across bases 2 through 5, and verify 864 scaled convolutions, 200 finite differences, and 2,100 lower triangular semigroup entries. Direct multiplicity enumeration independently checks 496 small dynamic-programming counts.

These finite checks establish agreement on the listed ranges. They do not prove the unbounded coefficient conjecture. The proposed settlement uses an induction on depth and argument to establish a scaled convolution identity, followed by matrix semigroup composition and row-recursive uniqueness.

The program is independently written for this report. The source rows and displayed cube are attributed to the OEIS entry.
