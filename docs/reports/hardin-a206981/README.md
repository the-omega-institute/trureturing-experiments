# Hardin A206981: exact finite checks

The target is the first FORMULA entry of R. H. Hardin's [OEIS A206981](https://oeis.org/A206981), revision 17 (19 June 2023), contributed 14 February 2012. Its NAME counts ordinary n by 2 binary arrays avoiding contiguous 010 and 101 in rows, columns, diagonals and antidiagonals. The offset is 1. The recurrence and generating function are marked empirical.

Run from the repository root:

```sh
python3 docs/reports/hardin-a206981/check.py
```

The recorded output is `check.out`. The program uses Python 3 integers only. It enumerates all arrays for heights 1 through 8 in all four literal directions, independently counts arrays through height 210 by a last-two-rows dynamic program, reproduces all 27 displayed OEIS terms, checks both switch-encoding inverses and the forbidden-pattern equivalence for every binary word of lengths 1 through 10, and checks the claimed recurrence and generating-series coefficients through height 210. It does not use the claimed recurrence to generate the comparison counts.

These are finite checks. They do not prove the all-length combinatorial equivalence, recurrence or generating function. There is no Lean proof in this report.

Program license: MIT. The sequence terms quoted in the program are from R. H. Hardin's OEIS A206981 and retain the [OEIS CC BY-SA 4.0 license](https://oeis.org/wiki/The_OEIS_End-User_License_Agreement). The remaining report text is MIT licensed.
