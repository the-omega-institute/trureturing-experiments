# Irwin's square-spiral hook-sum identity (OEIS A138099)

Source: [OEIS A138099](https://oeis.org/A138099), revision 26
(2022-09-08 08:45:33), sequence by Paul Curtz, conjecture by
Benedict W. J. Irwin (2016-05-14). The quoted triangle terms and three-line
diagram are short source excerpts; OEIS source text is governed by its
[End-User License Agreement](https://oeis.org/LICENSE).

Run with Python 3, using only the standard library:

```sh
python3 docs/reports/irwin-a138099/check.py
```

`check.out` is the output of that command. All calculations use integers.
The program reproduces all 77 published terms, the first eight displayed
triangle rows, Merca's formula through index 1001, and all ten displayed
spiral diagrams (including the empty diagram), up to translation. It builds
the spiral by unit steps west, south, east, north, with leg lengths
1,1,2,2,3,3,... . It counts ordinary hooks directly from cell coordinates,
then compares every growth increment through 1001 with the independently
flattened defining triangle. A separate pair count checks those totals;
all eight coordinate symmetries are also checked through 100 cells.

Every assertion must pass; the final line is `ALL_CHECKS_OK`.
These checks support source fidelity and finite consistency. They do not
prove the conjecture for every index, establish uniqueness of an implicit
source sequence, or certify the literature screen. The planned all-index
settlement uses a general pair-count and boundary-strip argument.
