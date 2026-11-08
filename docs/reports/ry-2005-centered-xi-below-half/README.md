# Centred Roynette–Yor xi transform below one half

The target is equation (V.35), p. 1275 of B. Roynette and M. Yor,
*Couples de Wald indéfiniment divisibles*, DOI 10.5802/aif.2125.
The proposed result excludes every nonnegative random variable satisfying
that centred transform identity when 0 < alpha < 1/2. The interval above
one half is outside the target.

Run `python3 check.py` with Python 3, SymPy and mpmath. The stored output is
`check.out`; successful execution ends in `ALL_CHECKS_OK`.

The exact checks reproduce the source's xi reflection identity at the
regular arguments -1 and 2, verify the conjugate-zero rational kernel and
a polynomial positivity certificate for its entire parameter range, and
verify the formal Taylor coefficient responsible for the negative second
difference. They also verify the square identity that every Laplace
transform must satisfy. Rational spot checks supplement those identities.

The mpmath section uses 100 decimal digits at alpha = 1/4. It prints
logarithmic derivatives and negative second differences at decreasing
positive t. These are numerical diagnostics, not exact or certified
bounds. The exact algebra does not certify the analytic zero sum,
its differentiation, the zero-height gap, or the Taylor remainder.
Those remain mathematical proof obligations for the proposed Lean work.

No proof of the target or its negation has been written in Lean here.
The proposed mechanism uses no RH assumption or moment hypothesis on the
putative random variable. The source and later distinction between
centred and uncentred transforms are discussed in Suvrit Sra,
arXiv:2611.02944v1, in the paragraph following the reciprocal-law proof.
