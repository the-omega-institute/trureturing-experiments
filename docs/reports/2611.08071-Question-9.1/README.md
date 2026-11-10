# Flag completion and normality

Original exact finite checks for Bruno Benedetti and Giuliamaria Menara,
*Diameter bounds for complexes that are close to flag or close to normal*,
arXiv:2611.08071v1, Question 9.1 (source submitted 5 November 2026).
Source: <https://arxiv.org/abs/2611.08071v1>.

Run with Python 3.12 or newer, using only its standard library:

```sh
python3 check.py
```

The recorded output is `output.txt`. All computations use finite sets, integers,
cardinalities, and exact graph searches. There is no floating point arithmetic.
The fidelity check reproduces the source's Section 8.3 example Gamma_7: its eight
listed facets give exactly the stated minimal nonfaces, it is normal, its dual
diameter is 3, and its completion is 4-dimensional.

For N = 7 and 8, the program constructs the boundary of the cyclic 4-polytope
from pairs of disjoint cycle edges, then relatively barycentrically subdivides
it while preserving its induced subcomplex on S = {0,1,3,5}. It constructs every
face directly from the relative-chain formula and independently compares the
result with descending stellar subdivisions of all non-A faces.

Every facet star is checked using the paper's induced-dual-graph definition,
including the empty face. Every link of codimension at least two is also
checked by its own vertices and edges. The program checks purity, every ridge's
degree, all cliques of the 1-skeleton, the exact completion identity, connectivity
of both full dual graphs, and the disconnected link of edge {3,5}. For d = 4 and
k = 2, the forbidden range 3 through d-k is empty.

The old core labels 0,1,3,5 are relabeled 0,1,2,3, respectively. Non-A original
faces are sorted by (cardinality, lexicographic vertex tuple), and their
auxiliary labels start at 4. Thus the printed bad edge (2,3) is the original
edge {3,5}. For N = 7 the f-vectors are (62,346,568,284) and
(62,346,570,285). For N = 8 they are (88,512,848,424) and
(88,512,850,425). In both cases the completion has exactly one bad facet star,
and the corresponding edge link has components of sizes 6 and 2.

These are finite experimental results. They do not certify the infinite-family
PL-sphere statement or supply a Lean proof. The intended formal refutation can
use the single N = 7 witness after preregistration. The implementation was
written from the mathematical definitions; no third-party program or Lean
source was copied or ported.
