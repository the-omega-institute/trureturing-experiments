# Hardin–Barker OEIS A209723: exact finite checks

Run `python3 docs/reports/hardin-barker-a209723/check.py` from the repository root. Python 3 and its standard library suffice. The stored output is `check-output.txt`; success ends with `ALL_CHECKS_OK`.

The program checks all 81 ternary squares and all 59,049 ordered pairs of five-entry rows using signed integer clockwise differences. It verifies the exact row graph as two three-vertex paths, eight isolated edges and 221 isolated vertices. An independent entry-by-entry enumeration checks heights n=1 through 8. Transfer counts divided by four match all 210 source b-file terms. It also checks both parity formulas, the recurrence and coefficients 0 through 210 of the cross-multiplied formal generating-function identity, using only integers and rational fractions.

These are finite fidelity and candidate checks. They are not an all-height proof or a Lean proof probe. The proposed all-height argument uses a source-array/walk bijection and the two-step continuation rule for a three-vertex path.

Source: [OEIS A209723](https://oeis.org/A209723/internal), revision 11 (15 July 2022); entry by R. H. Hardin, with Colin Barker’s FORMULA conjectures. The included `b209723.txt` is the official table from [Hardin’s b-file](https://oeis.org/A209723/b209723.txt), attributed to the OEIS Foundation and contributors under [CC BY-SA 4.0](https://creativecommons.org/licenses/by-sa/4.0/). Its SHA-256 is `67c1ee54560cfb5da53a7bdf7b21b4141856fff5321130ba8d1145d9d341305b`.
