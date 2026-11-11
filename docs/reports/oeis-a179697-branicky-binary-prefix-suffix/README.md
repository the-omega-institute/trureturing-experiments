# Binary prefix and suffix in OEIS A179697

The source is Dylan Hamilton's OEIS A179697, revision 23, with Michael S. Branicky's conjecture that every term starts with binary 1011 and ends with binary 11010. Binary reversal follows OEIS A030101, including deletion of leading zeros after reversal. Membership here explicitly means a positive integer satisfying both defining equations; zero satisfies the raw equations but is omitted by the published sequence.

Run with Python 3.12 or newer:

```sh
python3 docs/reports/oeis-a179697-branicky-binary-prefix-suffix/check.py
```

The standard-library-only program uses exact integer arithmetic. It reproduces A030101's example at 100 and initial values, and A179697's example at 90. It checks all 498 entries of the supplied table and exhaustively examines every integer from 1 through 10,000,000. It also checks the proposed normal form: k = 2p, where p is an odd binary palindrome and reverse(3p) = 5p. A second shift-based implementation checks reversal independently on the table and all inputs through 65,535. The saved output is check.out. These finite checks support a preregistration; they do not prove the universal conjecture.

Sources:

- [A179697](https://oeis.org/A179697), revision 23, 2026-11-08 14:55:24 (OEIS revision display); conjecture attributed to Michael S. Branicky, 2026-11-08; sequence attributed to Dylan Hamilton, 2010-07-24.
- [A030101](https://oeis.org/A030101), revision 203, 2026-06-28 15:47:50.
- [b179697.txt](https://oeis.org/A179697/b179697.txt), Michael S. Branicky, table n = 1..498. The adjacent b179697.txt is this source dataset, retained under the OEIS CC BY-SA 4.0 attribution and share-alike terms; see [OEIS license](https://oeis.org/LICENSE). It is not relicensed as program code.

The check program was independently written for this report and is Apache-2.0 licensed under the experiments repository's license.
