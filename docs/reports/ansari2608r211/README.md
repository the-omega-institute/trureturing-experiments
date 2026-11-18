# Ansari–Rockel upper inner curve: exact arithmetic check

The program reproduces the binary-model coefficient formulas in arXiv:2608.20176v1, Section 2.3, from exact integrals of its piecewise conditional CDFs. It independently integrates a joint rectangle density with uniform margins to obtain an unequal-cell counterexample family to upper-curve sharpness in Remark 2.11(a).

Run with Python 3.12 or newer and SymPy 1.14:

    python3 check.py

The stored output is `check.out`. All calculations use integers and exact symbolic rational arithmetic. The all-branch separating bound is checked through polynomial identities with explicit nonnegative factors on the stated domains. The interval certificate applies to every real parameter between 59/100 and 61/100. A positive-density variant is also checked. No sampled floating point calculation supports a conclusion.

The result is an arithmetic and analytic certificate, not a Lean proof. Constructing the measure and connecting its conditional statistics to the source definitions remain formalization obligations. The full exact region and the lower inner curve's sharpness are outside the settlement.
