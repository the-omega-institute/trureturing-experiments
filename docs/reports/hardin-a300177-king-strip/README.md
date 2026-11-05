# OEIS A300177: binary king-neighbour strips

`check.py` uses Python 3 integers exclusively. Run `python3 check.py`; its output is
stored in `output.txt`. OEIS A300177, revision 16, NAME, DATA and both empirical
FORMULA lines are the source: <https://oeis.org/A300177>. The parent entry is
<https://oeis.org/A300182> (column 3).

The program exhausts all anchored binary arrays of heights 1 through 6 using the
original cell-neighbour definition, and compares each array with avoidance of
three successive 000 rows or three successive 111 rows. It checks the three
terminal-run states against those literal counts, reproduces all 19 displayed
source terms, compares with a separate DP retaining all two-row symbol pairs
through height 210, and checks the matrix polynomial identity and generating
function coefficients using exact integer arithmetic.

These finite checks support the proposed uniform combinatorial argument. They
are neither a proof for all heights nor Lean kernel verification. The unbounded
argument requires an array-to-row-avoidance equivalence and append/delete
bijections for the three states.
