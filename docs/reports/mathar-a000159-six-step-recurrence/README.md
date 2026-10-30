# Mathar's A000159 recurrence

`check.py` checks the first FORMULA conjecture of R. J. Mathar in
[OEIS A000159, revision 54](https://oeis.org/A000159/internal).
The sequence is the coefficient of degree three of the ménage hit polynomial
in [OEIS A058087, revision 66](https://oeis.org/A058087/internal), attributed
there to Kaplansky and Riordan. The zero-row convention is the source's
explicit `U_0 = 1` convention. This report concerns the exact recurrence,
not van Hoeij's separate prime-shift congruence.

Run with Python 3 and SymPy:

```sh
python3 docs/reports/mathar-a000159-six-step-recurrence/check.py
```

The program uses only integers, rational numbers and exact symbolic algebra.
It reconstructs all six displayed initial hit polynomials and the first six
A000159 values from the finite rook sum. It checks the recurrence at every
index from 9 through 300, independently of recurrence-generated values.
It also compares a separate factorial-series expansion through degree 100
and checks the proposed rational differential certificate and its translation
to all seven recurrence coefficients. The saved output is `output.txt`.

Finite agreement alone does not prove the unbounded statement. The symbolic
checks verify algebraic identities used by the proposed formal-series
argument; no Lean proof or Lean probe is included. SymPy and the Python
runtime are outside the Lean trust boundary. The constants and initial
sequence data come from the cited OEIS entries; the checking program is an
independent implementation.
