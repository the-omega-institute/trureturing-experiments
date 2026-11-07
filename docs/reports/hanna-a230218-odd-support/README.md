# Hanna A230218: finite checks of the odd support

Paul D. Hanna's OEIS A230218 defines the normalized ordinary series by
`a(0) = a(1) = 1` and `[x^n] A(x)^(n^2+n+1) = 0` for `n > 1`.
Its first FORMULA conjectures that, for positive n, a(n) is odd exactly at
powers of two. The source is revision 18, with revision timestamp
2017-09-27 04:39:28: https://oeis.org/A230218/internal .

Run `python3 check.py` from this directory. The program requires only the
Python 3 standard library and uses integer arithmetic throughout. It independently
reconstructs the integer coefficients through degree 64, verifies their defining
rows, and reproduces the source's initial example. A separate binary polynomial
implementation reconstructs all reductions through degree 520 and compares them
with every archived OEIS term. It also checks the signed Catalan comparison-row
identity through degree 40, binomial parity through 300, and the two binary
valuation obstructions through 100000. Output is in `output.txt` and ends with
`ALL_CHECKS_OK`. These finite checks support a proposed all-index argument;
they do not constitute a proof or Lean verification of the conjecture.

The unmodified snapshot `b230218.txt` contains the 521 terms at
https://oeis.org/A230218/b230218.txt , attributed by OEIS to Paul D. Hanna.
Its SHA-256 is
`4f6facf9960662a653be9aa37c1a411d64128d405637d20bdfeecb7d100fdc89`.
The snapshot retains OEIS's Creative Commons Attribution Share Alike 4.0
licence: https://oeis.org/wiki/The_OEIS_End-User_License_Agreement .
The independently authored check program follows this repository's licence.
