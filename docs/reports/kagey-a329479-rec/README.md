# Kagey A329479 recurrence checks

`check.py` checks the degree-exactly-n interpretation of Peter Kagey's
[OEIS A329479](https://oeis.org/A329479/internal), revision 25, using Python
integers and SymPy exact rational algebra. It runs offline:

```sh
python3 docs/reports/kagey-a329479-rec/check.py
```

Requires Python 3 and SymPy. `check.out` was produced with Python 3.14.8 and
SymPy 1.14.0. Success means exit 0 and a final `ALL_CHECKS_OK` line.

Fidelity checks reproduce all 37 terms printed in the entry and exactly the
six degree-seven polynomials in its example. Exhaustive Horner evaluation
checks every monic binary polynomial at degrees 1 through 18. A separate
joint-evaluation dynamic program produces the counts through degree 1000
without assuming the conjectured recurrence; those counts agree with the
parity-binomial formula and satisfy the order-eight recurrence for every
index from 9 through 1000. Exact rational cancellation verifies the proposed
generating-function identity and its denominator.

The source's string-length comment has a one-position indexing mismatch with
its NAME and degree-seven example; the checks use the latter two. The
recurrence coefficients, initial values and example exponent sets are source
data. The check program is original research-seat code. No source fulltext or
third-party program is vendored here.

The finite checks and symbolic identities support a preregistration outline.
They do not provide a Lean proof of the polynomial-count bridge or of the
all-degree recurrence.
