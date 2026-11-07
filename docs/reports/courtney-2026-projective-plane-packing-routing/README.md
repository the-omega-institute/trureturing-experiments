# Whole-line packing routing in finite projective planes

The check uses Python 3 standard-library integer arithmetic. Run:

```sh
python3 check.py
```

The projective models PG(2,q) use normalized homogeneous coordinates over the prime fields q = 2, 3, 5, 7. The program verifies the incidence axioms, every point/line count, the source's weighted adjacency formula A = J - I and its exact eigenvectors, and the Ramanujan inequality by squaring nonnegative integers. These reproduce the adjacency/spectral identity in Section 2 and the Fano/PG(2,3) parameters of arXiv:2605.02498v3.

Full BFS on the Fano plane includes all 35 nonidentity permutations supported on a whole line. It visits all 5,040 permutations and has diameter 4, with layer counts 1, 35, 581, 2919, 1504.

A separate construction fixes a boundary line, constructs an affine incidence grid, factors the regular bipartite token multigraph into perfect matchings, and performs serialized row-column-row permutations. The check verifies every whole-line move, all Hall inequalities, preservation of the boundary and final token destinations. It checks all 5,040 Fano permutations and explicit deterministic test sets of sizes 1,092, 1,497 and 2,654 for q = 3, 5 and 7. Their maximum constructed costs are 9, 13, 21 and 29, respectively, within 4q+1 and hence 5q. Larger-order test sets are not exhaustive and do not test non-Desarguesian planes. The uniform claim needs the ordinary incidence/Hall argument and subsequent Lean proof.

PG(2,3) exact BFS is not run: its 6,227,020,800 states and 299 nonidentity generators entail 1,861,879,219,200 full expansion transitions; even a one-byte distance array takes 6,227,020,800 bytes and a full eight-byte queue takes 49,816,166,400 bytes, before representation overhead. These are exact counts for uncompressed full-state BFS, not claims excluding a compressed algorithm. No exact PG(2,3) routing number is reported.

The stored output is `output.txt`; success requires exit 0 and final line `ALL_CHECKS_OK`. No floating point, external package, numerical eigensolver or Lean proof is used.
