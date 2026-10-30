# Hanna A211898: exact finite checks

Primary source: [OEIS A211898](https://oeis.org/A211898/internal), revision 9,
with the exponent defined by [A089792](https://oeis.org/A089792/internal).
The source conjectures that the ordinary coefficients of

    exp(sum((2^k - (-1)^k)^k * x^k/k, k >= 1))

satisfy v_3(a(n)) = n - v_3(n!) for every n >= 0.

Run with Python 3, using only its standard library:

    python3 docs/reports/hanna-a211898/check.py

`output.txt` stores its exact output. The program reproduces all 13 coefficients
printed in the source, independently expands the exponential through degree 24,
and checks the coefficient recurrence, the factorial-scaled integer recurrence,
positivity, ordinary integrality, the exact valuation, and the involution residue
cycle through degree 120. The separate source formula a(n) = 3 modulo 6 is also
checked at indices 1 through 120. All evidence uses integers and rational numbers.

These finite checks test fidelity and support a proposed universal argument.
They do not prove the unbounded conjecture. No Lean proof is included.
