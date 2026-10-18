# A362534 exact table through n = 10,000

The full existing table of reduced fractions is in [bounds_1_10000.csv.gz](verification/bounds_1_10000.csv.gz). The table and [independent checker](verification/check.py) are byte-for-byte copies of the files in the 2026-10-19 verification archive. The archived public OEIS text export is included at [sources/entry-internal.txt](sources/entry-internal.txt), as required by the checker.

From this directory, run with Python 3 and its standard library:

```sh
shasum -a 256 -c SHA256SUMS.txt
python3 verification/check.py
```

The checker streams exactly 10,000 rows, checks every reduced-fraction field with exact arithmetic, and compares all 36 displayed OEIS terms. It uses an independent central-binomial recurrence and cross-checks selected values with `math.comb`. Run without `-O` or `-OO`; the checker requires assertions. It writes `verification/verification.json` and `verification/cancellation-factors.csv`.

CSV schema: `N,f_num,f_den,g_num,g_den,ratio_num,ratio_den,c`.

The finite check reports 396 exceptional indices in 198 odd/even pairs, starting at n = 19 with cancellation factor 5.
