# Rubin log-partition rigidity: exact finite checks

The check reproduces Proposition 3.1 of Michael P. Rubin, *A uniqueness theorem for exact separable log-partition decompositions*, arXiv:2607.22710v2. It also checks the canonical identity proposed for the bounded-log-ratio extension in Remark 7.4, including tightness, cancellation of a trial-density gauge, reverse tilts, and invariance under constant shifts of the log-weight.

Run with Python 3:

```sh
python3 docs/reports/rubin-log-partition-rigidity/check.py
```

The finite domain consists of every distinct strictly positive probability vector on 2, 3, or 4 atoms whose coordinates have a common denominator at most 5, and every exponential weight vector with coordinates in {1, 2, 3}. All arithmetic uses integers and `fractions.Fraction`. Logarithms of positive rationals are represented exactly as rational linear combinations of logarithms of primes. Nonnegativity of these expressions is checked by clearing coefficient denominators and comparing exact integer products. The program performs no floating-point calculations.

`check.output.txt` contains the program's output. There are 5,454 identity and gauge checks, 756 tightness checks, 206 reverse tilts, and 10,908 constant-shift checks. These checks do not prove uniqueness for unrestricted A and B, or the infinite-dimensional theorem. That theorem requires a general analytic argument; the finite checks only test fidelity and consequences of the proposed settlement.

The source is available at https://arxiv.org/abs/2607.22710v2. The source's assertion is not proved or refuted in Lean here.
