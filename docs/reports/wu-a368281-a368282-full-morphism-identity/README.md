# Wu's repunit-morphism partial-sum identity

`check.py` computes the fixed-word prefixes and numerical recurrence independently,
using only Python integers. Run `python3 check.py` from this directory with Python
3.12 or newer. The recorded result is `output.txt`; its final line is
`ALL_CHECKS_OK`. These finite checks do not prove the conjecture.

The program checks the identity and its two block-count cases for every denominator
1 through 64 and every index 0 through 10,000. It reproduces Table 1 of Jeffrey
Shallit, *Proof of Irvine's Conjecture via Mechanized Guessing*,
[arXiv:2310.14252v2](https://arxiv.org/abs/2310.14252v2), and compares against the
published OEIS data below. In particular S(10,000) for denominators 1 through 6 is
6180, 7320, 7914, 8284, 8541, 8730.

The target is Chai Wah Wu's FORMULA conjecture in
[A368281](https://oeis.org/A368281) and [A368282](https://oeis.org/A368282).
The first table indexes word positions from 1; the second indexes recurrence
values from 0. Both tables are read by downwards antidiagonals. All 87 entries of
each table b-file are checked. For row 6 these cover word positions 1 through 8
and recurrence indices 0 through 7. Separate recurrence b-files for rows 1
through 5 support checks through 10,000. No separate row-6 b-file was located;
row 6 beyond that table coverage is checked by the two independent constructions.

## Inputs and attribution

The files in `data/` are unmodified OEIS b-file snapshots retrieved 2026-11-11.
The program verifies each input's SHA-256 before using it. OEIS data remain under
[CC BY-SA 4.0](https://creativecommons.org/licenses/by-sa/4.0/), copyright
The OEIS Foundation Inc. and the individual contributors. Source URLs and offsets:

| File | Source | Meaning |
| --- | --- | --- |
| b005206.txt | https://oeis.org/A005206/b005206.txt | Hofstadter G, recurrence row 1, index 0 |
| b005614.txt | https://oeis.org/A005614/b005614.txt | Fibonacci word, row 1; file index 0 is source position 1 |
| b285431.txt | https://oeis.org/A285431/b285431.txt | Clark Kimberling's word, row 2, index 1 |
| b286389.txt | https://oeis.org/A286389/b286389.txt | Ilya Gutkovskiy's recurrence row 2, index 0 |
| b366860.txt | https://oeis.org/A366860/b366860.txt | Gutkovskiy's recurrence row 3, table by Amiram Eldar, index 0 |
| b366870.txt | https://oeis.org/A366870/b366870.txt | Gutkovskiy's recurrence row 4, table by Eldar, index 0 |
| b366871.txt | https://oeis.org/A366871/b366871.txt | Gutkovskiy's recurrence row 5, table by Eldar, index 0 |
| b368281.txt | https://oeis.org/A368281/b368281.txt | Wu's word table, flattened index 1 |
| b368282.txt | https://oeis.org/A368282/b368282.txt | Wu's recurrence table, flattened index 1 |

The independently written checking program is available under Apache-2.0.
