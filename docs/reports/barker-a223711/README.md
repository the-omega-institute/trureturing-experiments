# Exact checks for Hardin–Barker OEIS A223711

Run `python3 check.py` with Python 3. The program uses only the standard library,
integer arithmetic and `fractions.Fraction`. Its stored output is `check.out`.
Success requires exit 0 and the final line `ALL_CHECKS_OK`.

The source is [OEIS A223711](https://oeis.org/A223711), revision 9, attributed to
R. H. Hardin (2013), with Colin Barker's two FORMULA conjectures (2018). The
companion [A223716](https://oeis.org/A223716) supplies the same array definition.
The numerical source fixture `b223711.txt` is Hardin's table for indices 1–210,
retrieved from <https://oeis.org/A223711/b223711.txt>. The program also embeds the
28 terms displayed by revision 9. The source fixture contains numerical facts;
the sequence entry and its attribution remain at the linked OEIS pages.

The fidelity check reproduces all 28 displayed terms and all 210 b-file entries
from the printed formula and rational generating function, and checks the
printed empirical recurrence at indices 6–210. Independent literal enumeration
of all labeled binary n-by-2 arrays for n = 0–8 checks the source's row and column
conditions against a definition using an unmarked weak peak. These counts also
check the row-sum fiber sizes and the unique three-case block decomposition.
The empty array is auxiliary: its count is 1, while the source generating series
has constant coefficient 0 and starts at n = 1.

Exact polynomial cross-products check the proposed peak-class generating
function and partial fractions. These algebraic checks and the finite counts
support a proposed unbounded combinatorial argument. They do not prove the
all-n counting claim and are not Lean kernel certificates. The positive
settlement route must supply that general argument; bounded enumeration is
only a check.
