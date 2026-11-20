# Exponential sphere entropy and balanced stars

Run `python3 docs/reports/kraus13c51/check.py`. The recorded output is
`check.out`.

The check reproduces Proposition 4.1 on the seven-cycle in Kraus, Dehmer,
and Schutte (MATCH 70, 2013, pp. 885–900) and verifies that the tree in
Figure 5 minimizes entropy among the 19 balanced-arm types of order 20.
It compares exact distance sums with breadth-first traversal on small
balanced stars and brooms. At each of orders 1024, 2048, and 4096, it
checks the prescribed broom against every balanced-arm type with the
source’s exponential weights at k = 1.

All reported comparisons use enclosing rational intervals. The program
bounds exp(−1) by the alternating factorial series and bounds logarithms
by the atanh series after binary scaling. It rounds interval endpoints
outward on a fixed dyadic grid. The finite results check three orders of
the proposed family.

[Source paper](https://match.pmf.kg.ac.rs/electronic_versions/Match70/n3/match70n3_885-900.pdf)
