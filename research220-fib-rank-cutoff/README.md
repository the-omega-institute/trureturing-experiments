# Fibonacci rank mask below the block-exchange cutoff

This experiment combines the Fibonacci rank identity

\[
p\mid F_r \iff z(p)\mid r
\]

with the external necessary condition that a hypothetical least Robin
counterexample has exponent at least two at every prime
\(p\le 7\,608\,793\).  It enumerates primes \(p\) up to that cutoff,
computes the exact Fibonacci rank \(z(p)\), and retains rows for which
\(z(p)\) is a prime rank \(r\ge7\).

The output contains 25,370 distinct prime ranks.  The largest is
3,804,397, and every recorded rank has at least one witness prime below
the cutoff.  These rows give a finite conditional exclusion for the
affine candidates \(1+F_rg\): if \(p\mid F_r\), then
\(1+F_rg\not\equiv0\pmod p\), while the valuation condition would force
\(p^2\mid1+F_rg\).  The computation does not assert an infinite
exclusion, and it does not prove Robin or RH.

Build and run:

```sh
clang++ -O3 -std=c++17 rank_mask.cpp -o rank_mask
./rank_mask 7608793 > rank_mask.tsv
```

The first line is a summary; subsequent lines are `rank`, smallest witness
prime, and witness count.
