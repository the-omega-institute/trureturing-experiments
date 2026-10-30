# Cycle-star edge irregularity strength

Source: Umme Salma, H. M. Nagesh and Narahari N., *On edge irregularity
strength of cycle-star graphs*, Proyecciones Journal of Mathematics 43(3),
725–741 (2024), DOI [10.22199/issn.0717-6279-5801](https://doi.org/10.22199/issn.0717-6279-5801),
[arXiv:2405.12263v1](https://arxiv.org/abs/2405.12263v1), Conjecture 3.1.

Run with Python 3 and its standard library:

```sh
python3 docs/reports/salma-nagesh-narahari-cycle-star-irregularity/check.py
```

All arithmetic is integral. The graph uses exactly the source's vertex indices:
root 0, leaves 1 through m, and the cycle path m+1 through k+m−1.
Repeated vertex labels are allowed; all unordered edges must have distinct sums.

The program reproduces the explicit labeling in Theorem 3.2 for 40 leaf counts,
checks the source edge set against an independently indexed cycle with leaves,
and finds exact minimum label bounds by a separate complete backtracking solver
for 3 ≤ k ≤ 12 and 1 ≤ m ≤ 6. The solver uses degree and edge-count necessary
conditions to reject infeasible bounds and assigns the cycle before choosing
available leaf labels; it does not use the proposed construction formula.

The proposed formula is K = max(m+2, ceil((k+m+1)/2)). Its two explicit path
constructions are checked for all 4 ≤ k ≤ 160 and 1 ≤ m ≤ 80, including label
bounds, complete edge-weight distinctness, internal sum inventories and exact
reserve cardinalities. The source's piecewise formula is checked on the entire
corresponding range with k ≥ 8. The stored output is `output.txt` and ends with
`ALL_CHECKS_OK`.

These are finite exact checks, not a proof of the unbounded conjecture. No Lean
proof of the conjecture or of its negation is included.
