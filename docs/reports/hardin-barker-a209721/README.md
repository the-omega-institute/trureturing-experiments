# Hardin–Barker OEIS A209721: exact checks

The source is [OEIS A209721, revision 43](https://oeis.org/A209721/internal), whose NAME counts one quarter of the (n+1) by 3 arrays over {0,1,2} with pairwise distinct signed clockwise differences on every consecutive 2 by 2 square. Its FORMULA gives an empirical order-three recurrence and rational generating function.

Run `python3 docs/reports/hardin-barker-a209721/check.py` from the experiments repository root. Python 3 standard-library integers and `fractions.Fraction` supply all arithmetic. The program reads the accompanying source table and prints the stored JSON result followed by `ALL_CHECKS_OK`.

The program tests all 81 local squares and 729 ordered row pairs, checks the exact decomposition into two three-vertex paths, two two-vertex paths and 17 isolated rows, and computes unrestricted row walks through height parameter 210. All 210 source table terms agree, as do the parity formula, empirical recurrence and coefficients of the cross-multiplied generating function. Separate entry-by-entry enumeration covers all 729, 19683 and 531441 arrays at n=1,2,3. The raw counts begin 12,16,20,28,36,52,68,100.

These are finite exact checks, not a proof for unbounded n or a Lean verification. The proposed all-height argument uses an array/walk bijection and the fact that every vertex of a three-vertex path has two two-edge extensions.

`b209721.txt` is the source table by R. H. Hardin, downloaded from <https://oeis.org/A209721/b209721.txt>. OEIS content is attributed to the OEIS Foundation and its contributors under [CC BY-SA 4.0](https://oeis.org/wiki/The_OEIS_End-User_License_Agreement). `check.py` is an independently written check program. `check-output.txt` contains its exact output.
