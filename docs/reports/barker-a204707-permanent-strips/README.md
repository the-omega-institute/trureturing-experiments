# Barker A204707 permanent strips

Run `python3 check.py` with Python 3; only the standard library is needed. The stored result is `check.out`, whose final line is `ALL_CHECKS_OK`.

The checker uses integer permanents and determinants, rational arithmetic, and an exact pair representation of Q(sqrt(17)). It enumerates all labelled binary arrays for 1 <= n <= 5, verifies the source's seven examples, checks all 64 ordered row pairs and the two-class walk quotient, and compares exact transfer counts against every supplied source term. The rational generating-function identity and radical values are checked exactly. These are finite consistency checks; the all-height argument belongs to the proposed proof, and no Lean proof is supplied here.

Source inputs are R. H. Hardin's tables [A204707, n=1..210](https://oeis.org/A204707/b204707.txt) and [A183682, n=1..200](https://oeis.org/A183682/b183682.txt), downloaded for this check. The definitions and conjectures are [OEIS A204707](https://oeis.org/A204707/internal), revision 9, with Colin Barker's FORMULA conjectures, and [OEIS A183682](https://oeis.org/A183682/internal), revision 9. The embedded examples are transcribed from A204707. OEIS data attribution: R. H. Hardin, Colin Barker, and The OEIS Foundation Inc.; [OEIS license and attribution policy](https://oeis.org/wiki/The_OEIS_End-User_License_Agreement), CC BY-SA 4.0.

The source index starts at 1. The auxiliary active-row walk count at index 0 is 5; it is not the number of source arrays at index 0. Consecutive subblocks give a(1)=13, whereas allowing all column pairs gives 7. The prospective A204707 proof imports the separately preregistered A183682 result and adds the all-ones-array decomposition; this checker independently verifies their numerical agreement.
