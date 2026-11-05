# Han graph-substitution treewidth counterexample

Source: Yuhao Han, *Tensor Contraction with ZX Calculus, a Graph-Theoretic
Approach*, arXiv:2611.01890v1, Conjecture A.8; Definitions 3.2–3.3 and A.7.
Each loose edge is completed by its own fresh leaf. Widths and the subtraction
in A.8 are integer valued. The checker uses Python 3's standard library only.
Run `python3 check.py`; the stored output is `output.txt`.

The fidelity checks reproduce the source's complete-graph formula
(Example 3.8, orders 1 through 7) and nontrivial-tree value 1
(Example 3.5, paths and stars of orders 2 through 8). Exact treewidth is the
minimum over every elimination choice, with fill edges, using integer bitsets
and memoization. Each minimizing order is checked again with a set-based
elimination implementation. A separate checker verifies vertex coverage,
edge coverage, running intersection and width for the displayed path bags.

For every integer r ≥ 5, let A and B be disjoint r-cliques, with joining edges
only a1–b1 and a2–b2. Take H = G[A]. Replace H by the edge s1–s2, retaining
the attachments s1–b1 and s2–b2. The resulting graph G′ contains B and the
ear b1–s1–s2–b2. There are exactly two cut edges, with distinct endpoints on
both sides. H and S are connected, and deletion of any one vertex leaves
G and G′ connected.

The following path bags give upper bounds:

- G: A; {a1,a2,b1,b2}; B.
- Completed H: {a1,x1}; A; {a2,x2}.
- Completed S: {x1,s1}; {s1,s2}; {s2,x2}.
- G′: B; {b1,b2,s1}; {b2,s1,s2}.

Their widths are r−1, r−1, 1 and r−1 respectively. The matching lower bounds
come from the r-cliques A or B and the edge in S. A clique is contained in one
bag: its vertex-occurrence subtrees intersect pairwise; the Helly property of
subtrees of a tree gives a common bag. Explicitly root the tree, choose the
closest-to-root vertex in each occurrence subtree, and choose a deepest one x.
For any other subtree, a common vertex y with the selected subtree lies below
both chosen closest vertices. Their depths force x to lie between the other
closest vertex and y, and hence in that subtree. Thus x belongs to every
occurrence subtree. This supplies the required clique lower bound without
assuming that arbitrary pairwise-intersecting sets have common intersection.

Therefore (tw(G),tw(H),tw(S),tw(G′)) = (r−1,r−1,1,r−1). Conjecture A.8 would
require r−1 ≤ 3. The additive violation is r−4, unbounded as r grows.
The smallest refuting member is r=5. The preceding member r=4 has equality.

The program computes all four widths exactly for r=4 through 7, and checks
bags, clique certificates, substitution, independent leaves, distinct cut
endpoints and single-vertex deletion connectivity for r=4 through 40.
`INSTANCE_R5` contains every vertex, edge, bag and minimizing elimination order
for the smallest counterexample. Finite checks do not certify the unbounded
family. The argument above is mathematical exposition; no Lean proof or
kernel verification of the conjecture's negation is claimed.
