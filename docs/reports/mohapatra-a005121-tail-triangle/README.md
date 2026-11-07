# Mohapatra A005121 tail triangle: exact checks

Source: [OEIS A005121](https://oeis.org/A005121), revision 88, source timestamp
2026-11-05T17:38:02-05:00. The sequence name is “Number of ultradissimilarity
relations on an n-set.” The final FORMULA line attributes the tail-triangle
conjecture to Rajesh Kumar Mohapatra, Oct 30 2026. The source sequence is
computed here using its printed Lengyel recurrence with a(1) = 1.

Run with Python 3 and no additional packages:

```sh
python3 docs/reports/mohapatra-a005121-tail-triangle/check.py
```

`check.out` is the program's complete output. All computations use exact
integers. The checks reproduce all 19 published DATA terms and the source's
four strict partition chains at n = 3. Direct partition enumeration also checks
chain totals and Stirling counts through n = 6. An independent Stirling
inclusion-exclusion calculation checks 1 ≤ k ≤ n ≤ 12.

The program compares the printed triangle to OEIS through n = 12 (and all
19 displayed terms), compares its diagonal to the Lengyel recurrence through
n = 50, checks the row/Stirling identity through n = 50, and checks the finite
polynomial tail elimination for truncation degrees 1 through 12. The n = 1
initialization is handled separately from the diagonal sum rule.

These checks support indexing and the proposed proof mechanism. They do not
prove the all-index conjecture or its combinatorial interpretation in Lean.
The check program is original; no third-party program was copied or ported.
