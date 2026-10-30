# Oudra's A248183 square-root floor conjecture

Source: [OEIS A248183, revision 19](https://oeis.org/A248183/internal), NAME,
FORMULA and EXAMPLE. Clark Kimberling supplied the sequence; Ridouane Oudra
supplied the FORMULA conjecture. The source's PARI program starts its search
at the positive integer 1.

Run `python3 docs/reports/oudra-a248183/check.py`. The saved standard output is
`check.out`. Python 3's standard library suffices. All evidence uses integers,
`fractions.Fraction` and `math.isqrt`; no floating-point calculation is used.

The program reproduces all 68 displayed DATA terms and the exact example
`a(4) = 2` directly from the finite sum. It checks the proposed identity for
every integer from 2 through 10000, the telescoped remainder through 10000,
selected constants in the proposed stability interval, strict threshold
equalities and the interval's two endpoint obstructions. These are finite
checks supporting a preregistration, not a proof of the universal identity.

The published EXAMPLE row for index 5 prints 0.004 for 1/n²; the exact value
is 1/25 = 0.04. The check uses the exact definition and does not adopt that
decimal transcription error.

The source is attributed to its OEIS contributors; the displayed DATA terms
are used for source fidelity. OEIS licensing information is at
<https://oeis.org/wiki/The_OEIS_End-User_License_Agreement>.
