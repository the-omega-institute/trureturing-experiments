# Bala's A281428 eventual-divisibility conjecture

`check.py` checks the rational formal-series interpretation of Peter Bala's
conjecture 1 in [OEIS A281428](https://oeis.org/A281428), COMMENTS,
and Paul D. Hanna's even-index exponential generating function. The inspected
OEIS revision is 49, with timestamp April 14, 2026, 14:22:42; the conjecture is
attributed to Bala on March 3, 2026. The b-file is credited to Vaclav Kotesovec.
OEIS content is available under [CC BY-SA 4.0](https://oeis.org/wiki/The_OEIS_End-User_License_Agreement).

Run with Python 3, using only the standard library:

```sh
python3 check.py
```

To use an already downloaded exact b-file:

```sh
python3 check.py --b-file /path/to/b281428.txt
```

`output.txt` records the successful run. The b-file hash in the output identifies
the source bytes. The program independently constructs integer EGF coefficients
from the coupled polynomial ODE and rational ordinary coefficients by extracting
the constant-one square root of `sum binom(3n,n)*x^(2n)/3^n`. It compares all 201
source values, the displayed mod-385 example, the scalar recurrence, parity,
the source's algebraic equation and the literal cosh integral fixed point through
degree 400. Every ordinary coefficient checked has denominator supported on 2
and 3. For all 133 positive moduli at most 400 coprime to 6 it checks every index
from `ceil(N/2)` through 200: 13,400 divisibility comparisons.

All calculations use integers or `fractions.Fraction`; no floating point is used.
The checks are finite evidence. They do not prove the universal conjecture or
its negation, do not certify a minimal onset, and do not settle conjecture 2.
