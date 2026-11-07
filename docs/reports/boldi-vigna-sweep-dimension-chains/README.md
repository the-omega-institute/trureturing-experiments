# Sweep dimension of products of chains

Source: Paolo Boldi and Sebastiano Vigna, *The lattice of antichains of
intervals: Dyck paths, level sets, and sweep dimension*,
[arXiv:2610.41385v1](https://arxiv.org/abs/2610.41385v1), Definition 1,
Proposition 3, and Section 6, second closing bullet.

Run `python3 check.py`. The stored result is `output.txt`; success requires
exit code 0 and final line `ALL_CHECKS_OK`. Only Python 3's standard library
and integer arithmetic are used.

The fidelity checks exhaust linear orders and solve the resulting finite
set-cover problem, reproducing the source's sweep dimensions for chains
with 2 through 6 elements, antichains with 2 through 5 elements, and complete
bipartite posets of sizes (1,1), (1,2), (2,2), and (2,3). Singleton posets are
also checked directly from Definition 1, including its empty interval.

For products of nonempty chains with top ranks b_i and R = sum b_i, the
proposed answer is the sum of the rank sizes at floor(R/2) and
floor((R-1)/2) when R >= 1, and 1 when R = 0. The script constructs the
symmetric chain partition, checks every subset for order convexity, and
checks every interval-outsider separation by the constructed sweep orders.
It checks the matching two-rank obstruction using convex hulls. For at most
6 points it additionally exhausts all linear orders and computes the
minimum by exact set cover, independently of the construction.

The exhaustive interval checks cover nine displayed examples (up to 16
points), plus all other ordered tuples of 1 through 4 factor sizes in
1 through 4 whose product is at most 12: 127 additional cases. For Boolean
lattices with n = 0 through 12, the larger checks verify the partition,
saturation, symmetry, rank polynomial, singleton-sensitive sweep count,
and binomial identity; they do not enumerate all convex subsets when n > 4.

These are finite checks, not a Lean proof or an all-parameter computational
verification. The proposed proof uses general constructions and counting.
