# Distinct high-retention rows in finite Markov chains

The checker tests a proposed refutation of Saurav Jadhav, arXiv:2606.27073v1,
Section 10, Question 1. It does not contain a Lean proof.

Run with Python 3 and SymPy:

    python3 check.py

`check.py` uses only exact integers, fractions, symbolic polynomial algebra,
and rational enclosures of natural logarithms. The logarithm interval is
obtained from the positive atanh series with a proved geometric tail bound;
the displayed decimal endpoints are rounded outward using integer arithmetic.
No binary floating point is used as evidence.

Fidelity: the source's half-lazy star on five vertices has center retention
zero, leaf retention 1/3 and localization 1/2, checked with exact prime-log
coefficients. The proposed family is checked at m = 5, 6, 8, 10: rational
entries, row sums, symmetry, positivity, uniform stationarity, distinct rows,
maximizing block and total-variation separation. Rational polynomial bounds
are checked symbolically for every m >= 5 after substituting m = t + 5.

The eta intervals use the explicit uniform-on-A witness for the lower bound
and the analytic data-processing bound eta <= 1 for the upper bound. They
are not optimizer outputs, and the program does not independently prove data
processing. The infinite-family limits and the source's claimed implication
require the separate mathematical argument and subsequent kernel verification.

`output.txt` is the recorded successful run, ending in `ALL_CHECKS_OK`.
