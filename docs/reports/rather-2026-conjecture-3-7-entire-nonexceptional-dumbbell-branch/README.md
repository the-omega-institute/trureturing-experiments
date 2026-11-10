# Entire nonexceptional dumbbell NF orbit

Source: Bilal Ahmad Rather, *The NF-operator and the NF-Numbers of Simplicial
Complexes*, arXiv:2605.30781v1, Conjecture 3.7, its entire otherwise branch.

Run with Python 3.12 or newer:

```sh
python3 docs/reports/rather-2026-conjecture-3-7-entire-nonexceptional-dumbbell-branch/check.py
```

The original program represents **every** subset of the ambient vertices by an
integer bit mask. It computes facets by one-vertex extensions and computes NF by
avoiding every facet; the proposed closed-form recurrence does not generate the
computed orbit. All arithmetic and comparisons are exact integers and sets.

The output reproduces Theorems 3.3 and 3.5 for all 15 pairs
2 ≤ n ≤ m ≤ 6. For the 14 nonexceptional pairs it checks the D^5 entry, all 74
C_r facet classifications and all 74 descending/terminal transitions, the full
orbit through n+m+2, and minimality up to arbitrary permutation using triangle
and face-cardinality invariants. This includes every source-table parameter
pair up to interchange of the blocks. The largest ambient set has 12 vertices.

The separately reported B_(2,2)=P4 has first isomorphic return 1 and first
labelled return 2; an explicit permutation is checked. This known discrepancy
in the displayed exceptional branch is outside the proposed new contribution.

`check-output.txt` is the program's stdout. The final line is ALL_CHECKS_OK.
These finite checks corroborate a proposed universal argument; they do not
prove the unbounded conjecture and contain no Lean proof.
