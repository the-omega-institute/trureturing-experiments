# Hanna A395835–A395838: exact finite checks

Run `python3 check.py`; the successful output is `check.out`.

The program independently computes the formal equation by coefficient recursion
and by the exact Lagrange formula, reproduces every displayed coefficient for
n=2..6, and checks the four published b-file prefixes for n=2..30.
It also checks the rational-inverse derivative numerator as a polynomial identity
for these n and checks positivity on an exact rational grid. These are finite
checks supporting a preregistered analytic proof route, not a proof of infinitude.
No numerical approximation or proposed asymptotic cutoff is used as evidence.

The A395838 table includes the initial seed M(2)=3, then backward differences.
A395835 FORMULA (3) uses a forward-difference indexing inconsistent with that
table at n=2 (2 versus 3). Infinite length is unchanged by this finite seed and
index shift. No distribution or asymptotic conjecture is being checked.

Sources: Paul D. Hanna, [A395835](https://oeis.org/A395835) revision 48,
[A395836](https://oeis.org/A395836) revision 23,
[A395837](https://oeis.org/A395837) revision 9,
[A395838](https://oeis.org/A395838) revision 13;
`https://oeis.org/A<number>/b<number>.txt`, retrieved 2026-10-25.
The frozen prefixes in the program are attributed OEIS data, CC BY-SA 4.0:
https://oeis.org/wiki/The_OEIS_End-User_License_Agreement .
The checker code is supplied under the repository license; embedded OEIS data
remain under their stated source license.
