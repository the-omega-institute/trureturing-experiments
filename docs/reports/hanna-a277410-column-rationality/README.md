# Hanna A277410 column rationality: exact checks

Source: Paul D. Hanna, [OEIS A277410](https://oeis.org/A277410), revision 64,
COLUMN GENERATING FUNCTIONS conjecture of November 5, 2016; record timestamp
2017-01-06T07:02:08-05:00. The exponential normalization is the one printed in
EXAMPLE. The target uses columns indexed by the coefficient of y^k and ordinary
column variable z, including the zero coefficient at n = 0.

Run with Python 3 and SymPy 1.14.0:

```sh
python3 docs/reports/hanna-a277410-column-rationality/check.py
```

The checker independently iterates the source fixed-point equation in Q[y][t]
modulo t^13 and verifies the original compositional equation. It reproduces all
78 displayed triangle entries in rows 1 through 12. A separate exact
exponential-polynomial implementation constructs columns k = 0 through 6 using
the formal Taylor recurrence. It verifies their ODEs, zero-constant integrals,
initial vanishing, exponential degree bounds, denominator divisibility, and exact
numerator degrees. It reproduces every printed numerator P_0 through P_5.
Coefficient convolution checks each rational identity through 20 coefficients
beyond the predicted numerator degree.

All arithmetic supporting the checks is integer or rational arithmetic. The
stored output is in `output.txt`; success is exit 0 with the final line
`ALL_CHECKS_OK`. These are finite checks for preregistration, not a proof of the
uniform all-column conjecture or a Lean verification.
