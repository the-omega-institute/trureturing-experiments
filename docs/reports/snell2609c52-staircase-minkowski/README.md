# Snellman all-ones staircase Minkowski sum

The independent program checks Conjecture 52 and the following facet-count assertion in Jan Snellman, *Hilbert functions over the exterior algebra and f-vectors in higher rank I*, arXiv:2609.23605v1, submitted 20 September 2026. It uses only Python 3 standard-library integer and Fraction arithmetic.

Run `python3 docs/reports/snell2609c52-staircase-minkowski/check.py`. The stored output is `check.out`; successful execution ends in `ALL_CHECKS_OK`.

For every 2 ≤ n ≤ 7 and 1 ≤ d < n, the program constructs the source staircase vertices in dimension N=n+d and their pair sums. Independently, in the original x coordinates, it enumerates every choice of N−1 active inequalities from the conjectured H-description, adjoining x₁=1, and solves the resulting square system by rational Gaussian elimination. Every feasible solution lies among the source pair sums, and every pair sum satisfies the H-description. The exact finite-difference identities also checked by the program map that H-description into z≥0, Σz=2; hence it is bounded. These inclusions establish equality of the two finite-case polytopes. Affine ranks of their active faces give the facet counts. This uses the usual vertex characterization of a bounded finite half-space polyhedron, not a floating point hull computation.

Fidelity: the source’s Theorem 51 extremal matrix and vertex formula are reproduced throughout the grid. Example 53, n=5 and d=2, has affine dimension 6, 19 vertices and 9 facets. Its printed phrase “with the equation x₁=5 understood” conflicts with Definition 36 and Conjecture 52: every all-ones pair sum has x₁=1. The check flags this source typo explicitly and retains the conjecture’s x₁=1 convention.

The check examines 21 parameter pairs, 3,150 active systems, and 8,092 mass decompositions, including the zero-overlap-mass case. When n−d=1, the omitted coordinate inequality cuts a face of dimension N−3 rather than a facet. All cases pass. These bounded checks support the preregistration; they do not prove the assertion for arbitrary n,d and are not proposed as the positive Lean settlement.

The program was written independently from the mathematical definitions. No third-party Lean or program source was copied or adapted.
