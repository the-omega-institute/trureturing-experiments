# OEIS A390424: exact recurrence checks

Source: [OEIS A390424](https://oeis.org/A390424/internal), revision 33
(26 January 2026), with R. J. Mathar's conjectural recurrence. The defining
series is `g = 1 + x*g^3`, with constant coefficient one, and the target series
is `A = 1/(2-g^2)^2`.

Run `python3 check.py` with Python 3 and SymPy 1.14.0. The stored output is
[check.out](check.out). All arithmetic uses integers, rational numbers, or
symbolic polynomial identities.

The program constructs coefficients from the functional equation rather than
from the conjectured recurrence. It reproduces all 23 source coefficients,
checks the source's Pell/binomial formula through index 80, and checks the
recurrence at every index from 6 to 80. The residual at index 5 is -5952.

The program also verifies the derivative formulas and an exact second-order
inhomogeneous differential certificate, both by rational substitution and by
polynomial division modulo `x*t^3-t+1`. It checks the six coefficient
polynomials and the degree-five forcing polynomial. These exact identities
support an unbounded formal-series proof route. The finite coefficient checks
alone do not establish the unbounded recurrence, and no Lean proof is supplied.
