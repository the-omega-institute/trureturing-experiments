# OEIS A316553: exact finite checks

The independently written `check.py` uses Python 3.12 or newer and only integer
arithmetic from the standard library. Run `python3 check.py > check.out` in this
directory. Success is exit 0 and the final line `ALL_CHECKS_OK`.

The source is Andrew Howroyd, [OEIS A316553](https://oeis.org/A316553), revision 23
(2021-12-07), created 2018-07-06. It counts individual elements of exact order two
in SL(2,Z/nZ), starting at n=1. Its first Formula line conjectures
`a(n) = 2^(omega(n) + min(3, valuation(n, 2))) - 1`.
The table `oeis-b316553.txt` contains only the numeric rows of
[Antti Karttunen's b-file](https://oeis.org/A316553/b316553.txt), n=1..1023,
computed 2021-12-06. Source terms and quotations are attributed to their named
contributors and the OEIS Foundation under the
[OEIS End-User License Agreement (CC BY-SA 4.0)](https://oeis.org/LICENSE).
The Python checker is independently authored; it does not copy the source's
PARI or GAP programs. The numeric source table retains its source licence.

The checker verifies the source's displayed three matrices at n=2, the zero
count at n=1, and the equality -I=I modulo 2. It compares the full four-entry
matrix enumeration against the adjugate-pruned enumeration at every n=1..32,
and checks the equivalence between the original equations and the adjugate
conditions on every raw matrix at n=1..20. For every n=1..2000, it enumerates
all matrices of the shape (a,b;b',a) with 2b=2b'=0, tests the original determinant
and square, and compares the resulting count with the conjecture. These are
exhaustive searches of a mathematically justified necessary domain, not full
n^4 searches at every n. Neither the conjectured formula nor CRT is used to
construct the matrix sets. It also checks the scalar-root lifting description
at powers of two through 256 and 1511 instances of Israel's multiplicativity
statement. `check.out` gives the check results and all 2000 computed counts.

For completeness of the pruning: in SL2 the inverse is the adjugate. Thus
A²=I forces A=adj(A), giving d=a, 2b=2c=0. Conversely these conditions together
with a²-bc=1 give determinant one and A²=I by multiplication. This argument is
part of the mathematical interpretation of the program, not a Lean proof.

These finite checks support a preregistration. They do not establish the formula
at all positive moduli, and no Lean proof is included.
