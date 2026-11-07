# Hardin–Barker A185858: trace-avoiding binary arrays

Source: [OEIS A185858](https://oeis.org/A185858), revision 9. The source
counts `(n+2) × 3` binary arrays with unequal traces of neighboring contiguous
`3 × 3` subblocks, divided by 128. At width 3 there is one column of subblocks,
so only vertical neighbors occur. Traces are ordinary integer diagonal sums;
array symmetries are not identified.

Run with Python 3, using only its standard library:

```sh
python3 docs/reports/hardin-a185858-trace-arrays/check.py
```

`output.txt` contains the output. The program exhausts every binary array for
`1 ≤ n ≤ 5` and checks every trace-word fiber at these lengths. The valid
counts are `512, 2816, 15872, 88832, 498176`, giving source terms
`4, 22, 124, 694, 3892` after division by 128.

An independent row-append calculation retains the last two actual rows and
the previous trace. For `1 ≤ n ≤ 200`, it agrees with the weighted trace-word
count, has integral normalization, and satisfies the printed recurrence.
The first 20 terms agree with the source's displayed data. The generating
function product is checked through degree 200. The radical formula is
checked through index 200 using pairs of rational numbers in
`Q[s]/(s²−13)`, with no floating-point calculations.

These are finite exact checks, not an all-length proof or a Lean proof.
The target remains the count of actual arrays; a recurrence-generated
replacement sequence is never assumed to equal that count.
