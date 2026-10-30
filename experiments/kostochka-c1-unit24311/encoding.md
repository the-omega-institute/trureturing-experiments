# Exact finite correspondence

Source: Alexandr Kostochka et al., *Turán number of three disjoint cliques of
arbitrary sizes*, arXiv:2610.38430v1, Section 11 C1,
<https://arxiv.org/html/2610.38430v1#ThmConjecture1>, CC BY 4.0
<https://creativecommons.org/licenses/by/4.0/>. The source predicates are
paraphrased here. This retention revision clarifies the treatment of existing
isolated vertices and separates checked finite evidence from whole-source
claims.

Vertices are the actual integers 0,...,n-1. An edge is a distinct unordered
pair (u,v) with 0 <= u < v < n. Both endpoint bits are set, with no diagonal
bits. The independent replay additionally builds neighbor sets and checks
symmetry and irreflexivity directly. Its witness check uses vertex sets and
edge endpoints instead of the solver's adjacency-bit operations.

For sizes k, a packing consists of pairwise disjoint subsets of the respective
exact sizes, each internally edgeless. Cross-set edges are allowed. Partial
and full packings are separate existential quantifiers: the full packing may
replace every set. Permuting requested sizes together with their set labels
preserves existence. Sorting sizes to choose a recursion order imposes no
source restriction on earlier k_i.

The literal checker enumerates every independent subset of each required size
from the currently available vertices and recurses after removing that subset.
Its size order is fixed for the call. Every full packing appears in this
recursion after permuting its labels to that order. Only completed exhaustion
returns None; a NodeLimit exception means unknown. Retention runs this solver
only on the pre-existing small controls, not on the 30,000 saved inputs.

The existing capacity-coloring algorithm gives each packed class capacity k_i
and an unrestricted slack class capacity n - sum(k). Ordinary classes must be
independent; edges within the slack class are permitted. A packing induces such
an assignment by putting unused vertices in slack, and a completed assignment
gives a packing. The algorithm first assigns the existing nonisolated vertices.
It fills remaining capacities with the existing isolated vertices of the same
graph. It adds no graph vertices. Total remaining capacity equals the number
of those isolates; placing them introduces no edge obstruction.

For clique component orders c_j, a packed set uses at most one vertex from each
clique. A binary allocation matrix with row sums k_i and column sums <= c_j
therefore gives a packing by choosing distinct vertices for the rows in each
clique. The integral capacitated bipartite-flow criterion is

    sum(k_i for i in R) <= sum_j min(c_j, |R|)

for every row subset R. For a fixed row count r, its largest left side is the
sum of the r largest sizes, giving the saved prefix test. The independent
replay also checks all row-subset inequalities on the small controls and the
644 saved obstruction records. This criterion concerns disjoint unions of
cliques. No reduction of arbitrary source-admissible graphs to that family is
provided.

For the equality control, ell=4, k=(2,2,2,2), n=10, and G is K7 plus three
isolates. Its 21 edges equal (2ell-1)(n-sum(k)+1). Three disjoint independent
pairs exist. Four pairs could use at most four clique vertices and three
isolates, fewer than their required eight vertices. Thus this is an obstruction
at equality, outside C1's strict edge hypothesis. The separate fixed-partial
control uses one edge {0,1} on four vertices. The independent pair {2,3} cannot
be extended by the remaining edge, while choosing {0,2} and {1,3} gives the
required full packing.

The saved input generator uses integer square root. Reconstruction of the
fixed 30,000 records matches the previously reported complete stream digest;
the retention checker reads the materialized input stream directly. Neither
the generated parameter domain nor its finite positive witnesses replace the
whole source event. The ell-boundary guard, strict inequality, earlier-size
freedom, and full all-set rechoice remain in force.
