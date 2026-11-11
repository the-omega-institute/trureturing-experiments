# Hardin–Barker OEIS A230509: successor neighbours modulo eight

Source: [OEIS A230509](https://oeis.org/A230509), R. H. Hardin's entry
(21 October 2013), revision 12 (19 December 2025); Colin Barker's empirical
generating function (17 March 2018). The embedded table is Hardin's
[210-term table](https://oeis.org/A230509/b230509.txt). Numeric source data
are attributed to OEIS and R. H. Hardin. This program is independently written
and uses only the Python standard library.

Run `python3 docs/reports/hardin-a230509-mod8-strip/check.py` from the repository
root. `check.out` is its exact output. The final line is `ALL_CHECKS_OK`.
All arithmetic is integer arithmetic, including polynomial determinants.

The source counts n by 2 arrays of colours 0 through 7, with upper left colour
zero. Every cell needs a horizontal or vertical neighbour whose colour is its
own colour plus one modulo eight. Boundaries do not wrap.

The check independently enumerates all arrays through height four, reproducing
the two displayed source examples. A separate dynamic program uses actual
colours and two flags for whether the current bottom cells already have their
required neighbours. It matches every source term from height 1 through 210.

Subtracting the bottom left colour yields 22 reachable states (d,p,q).
Six states with both flags false have no accepting continuation. The live
classes are A={(1,1,0),(7,0,1)}, B={(3,1,0),(5,0,1)},
C={(3,0,1),(5,1,0)}, D={(1,1,1),(7,1,1)}, and
E={(d,p,q): d even and exactly one flag is true}.
The program checks transitions for every state and every absolute left colour.
With row vectors, initial vector (2,0,0,0,0), and terminal vector (0,0,0,1,0),
the live transition matrix is

    2  1  0  0  0
    2  0  1  0  0
    1  0  0  1  0
   14  2  2  2  8
    2  0  0  0  1

The exact determinant of I−xT is 1−5x+4x²+5x³−9x⁴−12x⁵.
Its adjugate entry (1,4) is x³(1−x), yielding the proposed numerator
2x⁴(1−x). The matrix annihilator and a signed integral conjugacy to an
all-ones terminal vector are also checked. The latter permits direct use of
the repository's generic scalarRelationShift theorem in a future proof.

The recurrence is checked at every height 6 through 210, and the generating
function identity coefficientwise through degree 210. These are exact finite
checks and an algebraic certificate for the proposed transfer. A proof for
all heights still requires the literal-array/transfer counting bridge and
induction. No Lean proof of the claim or its negation is included.
