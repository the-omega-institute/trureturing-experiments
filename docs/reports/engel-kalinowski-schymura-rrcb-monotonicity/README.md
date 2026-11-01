# Refined rectangle covering and LEDM restrictions

`check.py` checks the definitions and the two LEDM branches of the unnumbered Section 6 Future Work question in Engel, Kalinowski and Schymura, *Computing Lower Bounds on the Non-negative Rank via SAT Solvers*, [arXiv:2610.11950v1](https://arxiv.org/abs/2610.11950v1). It uses only Python 3's standard library and exact integers and finite sets.

Run from this directory:

```sh
python3 check.py
```

The stored result is `check.out`. Source Example 1 gives RCB 1 and RRCB 2, using two coincident rectangles with different labels. Source Example 2 gives RCB 2 and RRCB 3. Exhaustive minimum searches give LEDM values RCB = (0, 2, 3, 4, 4) and RRCB = (0, 2, 3, 4, 5) for sizes 1 through 5. All pairs of these sizes satisfy both monotonicity inequalities.

The search enumerates maximal positive-support rectangles and multisets of rectangle labels. Enlarging a rectangle preserves coverage and every touching-label requirement, so maximal rectangles suffice. Accepted covers are checked independently with finite sets. For all 81 matrices of size 2 by 2 with entries in {0, 1, 2}, minima computed using maximal rectangles agree with minima computed using all nonempty rectangles. Every refined family with at most two labels in that same matrix collection, including empty and repeated rectangles, is checked on every row/column restriction. The checks also include restrictions of the minimum LEDM covers, zero matrices, and a distance matrix with repeated coordinates. In total 12,664 family/restriction checks preserve the exact sets of touching labels.

These finite checks support source fidelity and expose errors in multiplicity or boundary conventions. They do not prove monotonicity for unbounded sizes. The planned mathematical argument restricts an arbitrary labelled cover to a submatrix and then compares least feasible sizes. Ordinary RCB monotonicity is prior art; the candidate new settlement concerns RRCB on the entire LEDM family. The regular-polygon branches are outside these checks and the proposed settlement.
