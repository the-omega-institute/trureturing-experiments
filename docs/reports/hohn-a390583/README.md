# Hohn A390583: exact repeated-radius checks

Run `python3 check.py` with Python 3. Only standard-library integer and rational arithmetic is used. The recorded stdout is in `output.txt`; success ends with `ALL_CHECKS_OK`.

The input `b390583.txt` is the 4000-term table from <https://oeis.org/A390583/b390583.txt>. The problem source is Charles L. Hohn, OEIS A390583, revision 10, published 2025, <https://oeis.org/search?q=id:A390583&fmt=text>. OEIS content terms: <https://oeis.org/LICENSE>.

The program reproduces the entry's first ten rows and checks the defining Descartes equation, ordering and common gcd for all 1000 data rows. It tests the conjectured cube restriction and the proposed primitive normalization on all 498 rows with repeated largest radius. The exact cancellation factors 1, 4 and 8 occur 278, 137 and 83 times.

An additional search covers every primitive repeated-largest row with positive first two radii at most 500, using the cleared equation to obtain the unique possible repeated radius. It finds 969 rows, with no failure. An independent gcd check covers all 76115 coprime pairs 1 ≤ u < v ≤ 500, including pairs outside the largest-radius condition. These finite checks do not prove the universal conjecture.

The entry's auxiliary family involving A000447 is erroneous as printed: at k = 3 it produces (3,5,420,420), which fails the defining equation, while the table gives (3,5,120,120). This auxiliary family is not used in the check or proposed proof. The program checks and reports this discrepancy explicitly. The nonprimitive row (3,6,72,72) provides a negative control for the common-gcd hypothesis.

Input SHA256: `0a0109ef4faf2b478342d8d38af2ed5bb6fa9744183f87d24553a8a147377084`.
