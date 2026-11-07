# Hermite Wronskian Mehler top coefficient

This report checks Conjecture 2 of Kuznetsov and Yuan, “Mehler formula for Wronskians of Hermite polynomials”, arXiv:2606.02405v1 (1 June 2026).

`check.py` uses SymPy’s exact integer and rational polynomial arithmetic. It implements the physicists’ Hermite polynomials and the ordinary derivative Wronskian. The source’s Proposition 1(iii), equation (7), recursively defines the coefficients `Q`; the checker verifies the coefficient identity from equations (1)–(2), reproduces Proposition 3 for `l=1`, checks the finite coefficient formula, and tests the Conjecture 2 top coefficient for every strictly increasing tuple of lengths 2, 3, and 4 drawn from `{1,...,6}` (50 tuples, including tuples beginning with 1).

The generated output is `check.txt`. These finite exact checks provide fidelity and small-case evidence only; they are not a universal proof.
