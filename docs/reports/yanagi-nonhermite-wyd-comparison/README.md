# Yanagi non-Hermite Wigner–Yanase–Dyson comparison

`check.py` uses only exact `fractions.Fraction` arithmetic. It evaluates the diagonal spectral formulas for `NV`, `NI`, and `NJ`, reproduces the explicit two-dimensional variance formula from Remark 2.1, checks `NJ = 2 NV - NI`, verifies the source bound `NI ≤ D/(M+m) NV`, and checks the factored squared gap for dimensions 2, 3, and 4. The exponent `alpha = 1/2` is used with rational square spectra, together with the endpoint cases `alpha = 0, 1`; the matrix entries are rational complex pairs. The final family `A = E_(1,n)` checks equality.

Run with Python 3:

```text
python3 check.py
```

The expected final line is `ALL_CHECKS_OK`.
