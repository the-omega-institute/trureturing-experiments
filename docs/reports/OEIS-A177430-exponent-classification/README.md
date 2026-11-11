# OEIS A177430 exponent classification

`check.py` reconstructs the least increasing square sequence from the coefficient recurrence

\[
n b_n=\sum_{j=1}^{n} a_j b_{n-j},\qquad b_0=1,
\]

for \(\exp(\sum_{n\ge1}a_nx^n/n)\). It reproduces the source values
\(a_1,\ldots,a_7=(1,9,16,25,36,144,169)\) and the exponential coefficients
\(1,1,5,10,24,51,122,244\) through degree seven.

It then evaluates the exponent family with exact `fractions.Fraction` arithmetic through degree 50 for every integer \(|m|\le20\), and checks that integrality occurs exactly for the nonnegative even exponents. It also checks that the first four coefficients have the same classification on this range. The final output line is `ALL_CHECKS_OK`.
