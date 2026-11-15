# Modified BIBD geodesic multiplicities

Source: Carlos E. Frasser, *Classifying K-Geodetic Graphs Via Combinatorial
Designs and Distance-regularity*, [arXiv:2611.05762v1](https://arxiv.org/abs/2611.05762v1),
Sections 6.9, 6.10 and 6.32, Problem 2.

Run with Python 3.12 or newer:

```sh
python3 docs/reports/frasser2611p2/check.py
```

`check.out` is the complete output. The program uses integers only. It builds the
point/flag graph from the indexed blocks and independently computes all distances
and shortest-path counts by breadth-first search. It compares those results with
the proposed complete table at every ordered vertex pair. The source fidelity
check reproduces Proposition 6.10.1: vertex number and the degrees of both types.

The named cases include the Fano 2-(7,3,1) design, a simple 2-(6,3,2) design,
the doubled indexed Fano design, all triples on six points, all pairs on four
points, and the boundary k=v with repeated indexed blocks. All simple designs
on 2 through 5 labelled points are checked, as are all designs on 2 through 4
points with each possible block occurring at most twice. Permutations of the
block indices give isomorphic graphs and are represented once in the exhaustive
part. The named and exhaustive sets overlap; 15,284 counts ordered pairs with
that overlap retained.

All eleven branches of the proposed table are exercised. The final line is
`ALL_CHECKS_OK`. These finite checks do not prove the general statement or
establish that the result is unpublished. Repeated blocks and k=v test an
extension whose exact inclusion in the source's terminology requires care:
the source calls the blocks a collection and describes incomplete as k<v.
