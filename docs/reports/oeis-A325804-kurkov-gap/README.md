# Kurkov's A325804 support gaps

The Python 3 standard-library checker uses exact integer arithmetic. It reproduces
the 58 displayed support positions in OEIS A325804, the 20-note table in
Drexler-Lemire--Shallit, and 35 values of A094373. It checks every adjacent
support gap for integers below 2^20, compares a separate bit fold and literal
product for integers below 4096, and exhaustively tests the proposed extremal
suffix formulas for 160 state/length pairs. These finite checks support a
proposed general proof; they do not prove the conjecture.

Run without external input:

```sh
python3 docs/reports/oeis-A325804-kurkov-gap/check.py
```

The stored output additionally compares every indexed row of the source b-file:

```sh
curl -L https://oeis.org/A325804/b325804.txt -o /tmp/b325804.txt
python3 docs/reports/oeis-A325804-kurkov-gap/check.py --bfile /tmp/b325804.txt
```

The b-file has 13,496 rows and ends at index 13,496 with value 65,535. Its link
text advertises 13,495 rows. The output records the downloaded file's SHA-256.

Sources: [A325804 revision 45](https://oeis.org/A325804/internal), including its
sole FORMULA conjecture and PARI convention P(0)=1;
[A004718](https://oeis.org/A004718/internal);
[A094373](https://oeis.org/A094373/internal); and
[arXiv:1402.3091v1](https://arxiv.org/abs/1402.3091v1), Introduction.
OEIS tables are credited to their entry authors under the
[OEIS license](https://oeis.org/wiki/The_OEIS_End-User_License_Agreement).
The check program is an independent implementation; the paper's text and code
are not redistributed here.
