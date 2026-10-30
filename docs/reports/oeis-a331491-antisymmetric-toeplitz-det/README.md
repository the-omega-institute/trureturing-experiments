# Spezia's antisymmetric Toeplitz determinant

Source: [OEIS A331491](https://oeis.org/A331491), revision 10, NAME, COMMENTS,
displayed example, PARI program and permanent data. The target is its determinant
conjecture, not an evaluation of the permanent sequence.

Run with Python 3; no third-party package is required:

```sh
python3 docs/reports/oeis-a331491-antisymmetric-toeplitz-det/check.py
```

`check.py` uses integers and exact fractions. It reproduces the displayed M(4),
its permanent 8, and the permanent data for n=0 through 4. It compares the literal
entry formula with an independent translation of the source's PARI program,
checks determinants for n=0 through 32 using Bareiss elimination, and independently
cross-checks n=0 through 8 using rational Gaussian elimination. For n=1 through 16
it checks every entry of the difference/sign reduction, the alternating null
vector, the border isolation, and the determinant-one splitting of the constant
skew matrix. All divisions in Bareiss elimination are checked to be exact.

`output.txt` contains the successful output, ending in `ALL_CHECKS_OK`.
These finite checks support a preregistered proof outline; they do not prove
the identity at all sizes and are not Lean verification.
