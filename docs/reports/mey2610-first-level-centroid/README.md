# First-level centroid covering

The question is the paragraph after Theorem 1.2 of Mark Meyer, *Finite
convexification by Minkowski summation*, arXiv:2610.26744v1. It asks whether
uniform boundedness of the first centroid-centered exposed-set covering radii
suffices without the polytope hypotheses in Theorem 1.3.

Run `python3 docs/reports/mey2610-first-level-centroid/check.py` from the
repository root. The program uses only Python's standard library and exact
integers and fractions. Its saved output is `output.txt` and ends with
`ALL_CHECKS_OK`. The original check program is licensed Apache-2.0.

The checks reproduce the exposed-interval calculation in the source's §2.3:
for a_j = 2^(-j-1) and b_j = 2^(-2j-1), the radius is a_j/b_j = 2^j.
They also check the source's uniform tail sums and 32 pairwise nonparallel
rational unit directions. This checks the infinite example's exact formulas;
it does not replace the source's infinite bodies with finite truncations or
certify compactness by finite computation.

The settlement checks cover the following explicit families:

- All 2,793 Cartesian products of nonempty subsets of {0,1,2} in dimensions
  2, 3 and 4, with rational axis-parallel box convexifiers. All 201,390 proper
  coordinate faces are examined; the product radius is independently checked
  by full nearest-point grid evaluation for the 343 three-dimensional pairs.
  The orthogonal-face radii satisfy both the source's polytope bound R_j≤R_1
  and the proposed bound R_j≤((n−1)!/(n−j)!)R_1.
- All 255 nonempty subsets of the three-dimensional unit cube's vertices,
  with a planar unit disc or vertical unit segment as convexifier. Rational
  directions certify infinite first-level radii for every nonsingleton set
  with the disc, and for all sets not on a common vertical line with the
  segment. The 12 remaining segment cases have R_2=R_1.
- In three dimensions, A = circle(a) × [-h,h] with K = disc(b) × {0}:
  R_1 = c = a/b and R_2 = 0. Also A = disc(a) × {-h,h} with vertical
  K = [-b,b]: R_1 = R_2 = c = h/b. All 125 triples from
  {1/3,1/2,1,2,3} are checked in each family. The geometric reductions to
  radial and vertical intervals are explicit analytic descriptions, and the
  program checks their exact thresholds; it does not enumerate the sets.
- For C = conv(unit disc in z=0, (1,1,0)), u=e_1, v=−e_2, x=(1,0,0),
  x is an endpoint of C_u that is not an exposed singleton of C. At 62
  rational parameters, the program verifies the projection construction's
  supporting-functional certificate and error ||z−x||²=4/(k²+1).
- The centroid difference-body inclusion G−G⊆(d+1)(G−g(G)) is checked at
  all 285 vertex pairs of standard simplices in dimensions 0 through 8.

These are finite tests of specified families and exact algebraic formulas.
They are neither a proof of the universal assertion nor a Lean verification.
The proposed general proof uses projection limits, compactness, the classical
centroid asymmetry inequality, and the source's Theorem 1.2.
