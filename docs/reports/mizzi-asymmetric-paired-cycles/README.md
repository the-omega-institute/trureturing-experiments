# Asymmetric paired cycles in unstable graphs

The dependency-free Python 3 program checks exact finite consequences of the proposed structural argument for the second assertion of the unnumbered Conjecture in Section 7 of Russell Mizzi, *Lifting and Folding: A Framework for Unstable Graphs and TF-Cousins*, arXiv:2603.27559v3. These results do not prove the universal assertion.

Run `python3 check.py`; the stored standard output is `output.txt`.

The fidelity checks reconstruct the twelve-vertex graph directly from Figure `gammacockroach` in the source TeX. They verify its stated TF pair, connectedness, distinct open neighborhoods, trivial ordinary automorphism group, and its paired cycle lengths for k = 3 and k = 5. A disjoint triangle and hexagon also occur. The k = 5 witnesses overlap, as fifteen distinct vertices cannot fit in this graph. The source's seven-vertex excluded example has exactly cycle lengths 5 and 6 and automorphism group order 4. Srivastava's arXiv:2608.15281v1 counterexample to the older unrestricted clause has the ordinary automorphism (1 6)(2 3), excluding it from the revised asymmetric claim.

The structural checks enumerate every graph invariant under the specified inverse TF pair for twelve prime-orbit layouts: p = 3, one to three moved orbits and zero to two fixed vertices, restricted to at most ten vertices; p = 5, one to two moved orbits and zero to one fixed vertex. The code generates all admissible edge-orbits, checks both adjacency and nonadjacency, independent moved orbits, uniform fixed-vertex adjacency, and the ordinary automorphism induced by each bipartite quotient. For every nonbipartite quotient it constructs the paired cycles and verifies their edges and injectivity. This covers 5,078 graphs. It does not enumerate all graphs of these orders or all TF actions.

The reflection check enumerates all matching-sum tuples for (p,k) = (3,3), (3,5), (3,7), (5,3), (5,5), and (7,3): 6,050 cases. For each it independently computes graph components and the return permutation and verifies one k-cycle and (p−1)/2 disjoint 2k-cycles. All calculations use integers and exact finite sets. The final line is `ALL_CHECKS_OK`.

Sources: https://arxiv.org/abs/2603.27559v3 and https://arxiv.org/abs/2608.15281v1. Source diagrams are used only to extract finite edge data; no TeX, diagrams, or article prose are redistributed here.
