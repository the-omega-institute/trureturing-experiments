# Kunisky spiked eigenvector fluctuations

`check.py` uses only integers and `fractions.Fraction`.  On the infinite
subsequence `N=m^4`, it verifies the exact isotropy identities

\[
\mathbb E W_N^2=I_{N+2},\qquad
\|\operatorname{Cov}(\operatorname{vec}W_N)\|=2N^{-1/4},
\]

for the proposed counterfamily, including positivity of all variance-profile
parameters and the unit norm of the nominated observation row.  It then checks
the exact limiting arithmetic from the optimizer and eigenvector calculations:
the predicted covariance is `1`, while the actual fluctuation variance is
`7/16`.  The script is a deterministic arithmetic check; it does not sample
the Gaussian matrices or claim to prove the random-matrix limit.

Run `python3 check.py`.  The final line is `ALL_CHECKS_OK`.
