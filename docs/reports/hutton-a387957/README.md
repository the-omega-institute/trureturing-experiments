# Hutton A387957: exact finite checks

Run `python3 check.py` with Python 3. The program uses the standard library and integers only. Its output is in `check.out`.

The source is Thomas Hutton's [OEIS A387957](https://oeis.org/A387957/internal), revision 20 (17 November 2025), whose two FORMULA lines call the generating function conjectured. Equal parts with equal colors are indistinguishable, and the color sum is 1 modulo 4.

The checks reproduce the three explicitly listed partitions at n = 2 and all 57 values in [Sean A. Irvine's b-file](https://oeis.org/A387957/b387957.txt), n = 0 through 56. Attribution for the quoted numerical table: OEIS Foundation Inc. and Sean A. Irvine; OEIS data are available under [CC BY-SA 4.0](https://oeis.org/LICENSE).

An unbounded-multiplicity dynamic program counts all four color-sum residues. A separate Euler pentagonal recurrence and ordinary convolutions compute the two- and four-colored partition counts. Literal multiset enumeration cross-checks both algorithms for n = 0 through 8. For every n = 0 through 128, the program verifies color reversal, the unsigned and signed counts, and the conjectured identity

    4*a(n) = p4(n) - (p2(n/2) if n is even, else 0).

The final output line is `ALL_CHECKS_OK`. These finite checks support source fidelity and the proposed route; they do not establish the assertion for every natural number. The planned proof uses a color-reversal bijection and a general coefficient argument, rather than positive bounded enumeration.
