# Spezia A001752 exact checks

Source: [OEIS A001752](https://oeis.org/A001752), COMMENTS conjecture by Stefano Spezia,
January 12, 2022; entry revision 119, August 30, 2026, offset 0.
The source defines the sequence by `1/((1+x)*(1-x)^5)`.

Run with Python 3 and SymPy:

```sh
python3 check.py > output.txt
```

The program reproduces all 43 published terms using formal reciprocal coefficient
recursion, and compares the source's floor formula and alternating binomial formula
through index 200. It computes characteristic polynomials over the integers for
matrix sizes 2 through 32, checks each principal minor of size 2, and compares the
second coefficient with the distance sum and the sequence. It separately checks
the distance sum through size 203. Python integer division supplies mathematical
floor for negative matrix entries. All evidence uses exact arithmetic.

These finite checks support the proposed settlement; they do not prove the
universal conjecture. No Lean proof is supplied.
