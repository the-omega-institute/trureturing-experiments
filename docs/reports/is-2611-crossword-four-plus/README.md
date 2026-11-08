# Four-plus crossword clue counts

Exact checks supporting a proposed proof of Isham and Schnitzer,
*Maximizing Crossword Grids with a Fixed Minimum Clue Length*,
arXiv:2611.03364v1, Conjecture 1.4.

Run `python3 check.py`. Python's standard library suffices. Every evidence
calculation uses integers or `fractions.Fraction`; no floating point is used.
The stored output is `output.txt`.

The program reproduces the displayed nine-by-nine grid's 26 clues, the
all-white cases four through eight, and the source's known values 80, 144,
228, 332 at sizes 15, 20, 25, 30. It constructs the proposed seed-and-extension
family for every multiple of five from 15 through 150 and verifies its
clue count, minimum run length, rotational symmetry and white-cell
connectivity directly.

At size 15, it enumerates all 32768 binary lines, the 296 lines satisfying
Four+, and all 87616 ordered pairs of legal lines. It exhaustively checks
the boundary-corrected pair inequality and minimizes a relaxation that
contains every Four+ grid. The relaxation proves the upper bound 80:
row-T plus column-S deficits are at least five, and transposition gives
another five. The program does not enumerate all grids or assume they
are connected or symmetric. The construction attains 80.

The check also enumerates the proposed summation's integer cases through
k=100. These finite checks support the proposed argument; they do not
constitute a proof of the unbounded conjecture. No claim or negation is
proved in Lean here.

Sources: <https://arxiv.org/abs/2611.03364v1> and
<https://oeis.org/A378174>. The check program is independently written
from mathematical definitions; no third-party program source is used.
