# Hanna square-exponent integrality check

`check.py` uses only integer arithmetic and `fractions.Fraction`.  It solves the
coefficient recurrence for the three cutoffs `M(n)=n`, `n+1`, and `2*n`, checks
the first twelve coefficients against the OEIS values where they are listed,
and checks the last-numerator divisibility by `M(n)` for every `n=2..13`.
The A350412 normalization is handled as `B=A-1`, so the reported coefficients
start with `a(1)` as in OEIS.  `check.out` is the exact captured output; the
final line is `ALL_CHECKS_OK`.
