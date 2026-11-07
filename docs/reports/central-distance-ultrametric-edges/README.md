# Central-distance ultrametric edges

Exact-arithmetic checks supporting a preregistration draft for Dovgoshey–Rovenska,
*On the Graph Corresponding to Center of Distances of Metric Space*,
arXiv:2610.11596v1, Conjecture 4.2.

Run with Python 3 (standard library only):

```sh
python3 docs/reports/central-distance-ultrametric-edges/check.py
```

`output.txt` contains the output. Every distance is an integer or a
`fractions.Fraction`; no floating-point calculation supplies evidence.

The center is the intersection of full row-distance sets, including zero.
Edges are unordered pairs of distinct vertices whose distance belongs to the
center. The fidelity checks exhaust all symmetric tables with positive
entries in {1,2,3} on 2, 3 and 4 vertices. Among the metric tables, they
reproduce Theorem 3.10's lower bound and sharpness. Among the ultrametric
tables, they reproduce Theorem 3.13's bound and its equality iff star condition.

The explicit four-point pair has central edge counts 6 and 4. The first space
is ultrametric and satisfies Conjecture 4.2(i); the second is an unrestricted
metric with a three-element center. The family checks cover every n from 4
through 32: the proposed ultrametric has 3n−6 central edges, while the
comparator has 2 ceil(n/2). Each comparator uses two disjoint spanning edge
covers, with distinct larger rational distances on the remaining pairs.
The program checks every triangle, exact centers, edge counts, both cover
conditions, and the exhibited partition and stars.

These are finite experiments and a concrete counterexample check. The
unbounded family and sharp minimum require the general argument in the
preregistration draft; the program does not certify a Lean proof.
