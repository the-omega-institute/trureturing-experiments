# Kurkov's A102086 row-polynomial identity

Source: [OEIS A102086](https://oeis.org/A102086/internal), revision 19,
8 July 2026; the conjecture is attributed to Mikhail Kurkov, 3 July 2026.
The matrix is Paul D. Hanna's triangle. Indices start at zero.

Run `python3 docs/reports/kurkov-a102086/check.py` from the repository root.
The stored output is `check.out`. The program uses only Python integer arithmetic.
It separately implements the matrix-square and polynomial recurrences, reproduces
all 47 published terms and all five matrix-square example rows, and compares every
coefficient through row 30. It also checks the proposed first-return identity
through row 30 and enumerates individual weighted terminating paths through row 7.

The polynomial diagonal is exactly X^n; the zero-column boundary is 1, including
P(0,0)=1. No row reversal or alternate diagonal is used. These finite checks
support source fidelity and the proposed route; they do not prove the conjecture.
