# Quantitative balancing in higher-dimensional product grids

Source: Jiechen Zhang, *Extremal Spanning Trees in Product Grid Graphs*,
arXiv:2606.24016v2, Section 6, “Quantitative balancing in higher dimensions”.
The question asks for a positive dimension-dependent coefficient controlling
logarithmic spanning-tree gain by the drop of half the lattice edge boundary,
for all equal-volume positive integer tuples in logarithmic majorization order.

Run with Python 3 and SymPy:

```sh
python3 check.py
```

`check.out` is the output, produced using Python 3.14.8 and SymPy 1.14.0.
The check exits zero and ends in `ALL_CHECKS_OK`. Its evidence uses integer
and rational arithmetic throughout; even the printed decimal intervals are
outward-rounded rational enclosures, not floating-point estimates.

The source fidelity check reproduces the Section 6 periodic quotient count
12321 and the standard 3-by-3 torus count 11664. These are independent of the
proposed quantitative statement. The program counts spanning trees using the
determinant of a reduced combinatorial Laplacian. Rational Schur complements
compute the determinant exactly. SymPy Bareiss determinants independently
cross-check five boxes of volumes 1, 6, 8, 12 and 24.

For dimension three and every volume from 1 through 200, it enumerates all
778 sorted positive triples (3722 triples before coordinate permutations).
It checks every logarithmically comparable pair: 2234 strict sorted pairs,
representing 58761 ordered strict pairs, and 778 sorted equality cases.
Coordinate permutations give isomorphic graphs and leave the boundary and
majorization unchanged. Prefix products of descending integer side lengths
check log-majorization exactly.

With D = sigma(n) - sigma(m) > 0 and c_3 = log(5/4)/2, the proposed inequality
is equivalent to the integer comparison

```text
tau(m)^2 * 4^D >= tau(n)^2 * 5^D.
```

Every checked pair satisfies it. The minimum ratio log(tau(m)/tau(n))/D occurs
at n=(3,8,8), m=(4,6,8), volume 192 and D=8. Its certified rational enclosure,
printed in decimal, is [0.515234586178623, 0.515234586178624]. The proposed
c_3 lies in [0.111571775657104, 0.111571775657105]. Cross-powers of exact
integers determine the minimum before any logarithmic enclosure is computed.
The enclosure uses the atanh series with an explicit geometric tail bound
and reduction to arguments between one and two.

The all-one box and its empty Laplacian cofactor are included with tree count
one. Zero boundary deficits occur precisely for identical sorted tuples in
the checked range. The two nontrivial example counts supplied with the
candidate are also reproduced exactly.

This is finite experimental evidence, neither a proof of the universal
assertion nor a Lean kernel verification. The proposed general argument uses
a curvature estimate for the path heat-trace interpolation; enumeration is
only a diagnostic helper and is not the proposed settlement method.
