# Spezia Pell circulants

`check.py` uses Python 3's standard library and exact integer and rational arithmetic.
Run `python3 check.py`; the stored output is `check.out`.

The source is Stefano Spezia's OEIS A392131 revision 16 (determinant) and A392132
revision 10 (permanent). Their Mathematica programs start with a once-right-rotated
Pell vector. The example at n=4 has determinant 640 and permanent 772. The source
asserts that permanent/determinant tends to 1.

The program reproduces all published terms at indices 0..13 and computes both
matrix functions through n=16. It cross-checks the permanent by subset dynamic
programming, Ryser's formula (through n=9), and direct permutations (through n=8),
and the Bareiss determinant by direct signed permutations (through n=8).
It checks the proposed rational entry majorant through Pell index 200, enumerates
all simple directed cycles through n=8, and verifies weighted composition counts
through total 17. Tail estimates at five specified sizes use fractions only.

These are finite checks of source fidelity and the proposed proof mechanism.
They do not certify the infinite limit. No Lean proof is included.

Sources: https://oeis.org/A392131/internal and https://oeis.org/A392132/internal.
The embedded numerical terms are factual sequence data; this program implements
independent algorithms and does not copy the source's Mathematica code.
