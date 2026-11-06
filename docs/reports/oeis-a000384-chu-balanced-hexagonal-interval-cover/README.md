# OEIS A000384 balanced hexagonal interval cover

`check.py` uses only exact Python integers and finite sets.  For every `1 <= n <= 30`
it constructs `N = n(2n-1)`, `A = {0,...,n-1} ∪ {N-n+1,...,N}` and
`B = {0,n,...,(2n-1)n}`, then checks `A,B ⊆ {0,...,N}`, both cardinalities `2n`,
and `A+B = {0,...,2N}`.  It also exhaustively enumerates every equal-size pair
for `n = 1,2,3` and reports the minimum `k`, and checks the counting inequality
`(2n-1)^2 < 2N+1` through `n = 30`.

The stored output is the exact stdout from Python 3.  The finite checks support
the preregistration and do not replace the proposed all-`n` proof.
