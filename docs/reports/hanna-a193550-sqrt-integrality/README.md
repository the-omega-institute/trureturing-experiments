# Hanna A193550: exact square-root integrality checks

Source: Paul D. Hanna, [OEIS A193550](https://oeis.org/A193550/internal),
revision 19, last edited August 25, 2017; entry authored July 30, 2011.
The sole COMMENTS conjecture asks whether the normalized formal square root
of `(A(x) - 1)/x` is an integer series. The continued-fraction weights are
`beta(j) = (2*floor(j/2) + 1)^2`, starting with `1, 1, 9, 9, 25, 25`.

Run with Python 3, without third-party packages:

```sh
python3 docs/reports/hanna-a193550-sqrt-integrality/check.py
```

`output.txt` contains the successful output. All arithmetic is exact.
The program reproduces the source's 17 displayed A coefficients, its 10
explicit B coefficients, and its observation that the first 500 B coefficients
are integral. It computes A by weighted Dyck paths and independently checks
direct finite nested fractions through degree 40, including two depths and
a nonconstant normalized terminal tail. It checks the Catalan congruence
modulo 8 through degree 500 and the proposed explicit integral lift through
degree 40. The displayed A and B coefficients are literal source data.

The root recurrence checks every division by 2 for zero remainder. A failed
check raises an exception; success ends with `ALL_CHECKS_OK`. These are finite
checks, not a proof of the conjecture or of infinite stabilization. No Lean
proof is included.
