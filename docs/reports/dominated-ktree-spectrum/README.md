# Dominated chromatic spectrum of k-trees

Source: Sneha R, J. Geetha and K. Somasundaram, *Dominated Coloring of Some
Classes of Perfect and Unicyclic Graphs*, arXiv:2609.32439v1, Section 2.3,
Open Problems item 1. Source: <https://arxiv.org/abs/2609.32439v1>.

Run `python3 docs/reports/dominated-ktree-spectrum/check.py` from the repository
root. The program uses only the Python standard library and integer arithmetic.
Its complete output is in `output.txt`; successful execution ends with
`ALL_CHECKS_OK` and exit status 0.

A color class must be independent and have a common adjacent vertex in the
**open** neighborhood. Every used class has this requirement. The program
exhaustively enumerates all such subsets and computes the minimum partition by
an exact subset dynamic program. It reconstructs and validates an optimal
partition. A separate restricted-growth partition enumeration verifies all
22 cases with at most eight vertices, without sharing the subset tables or DP.

Fidelity checks reproduce the source's path formula for orders 2 through 16:
`floor(n/2)` when `n` is divisible by four, and `floor(n/2)+1` otherwise.
The source's maximum-degree lower bound is checked on every tested graph.

The realization family has vertices `0,...,n-1`, with distinct vertices adjacent
exactly when their distance is at most `k`. All 54 cases with `k=1,...,4` and
`k+1 <= n <= 16` are checked. Actual deletion steps verify k-tree membership;
connectivity, the initial clique value, degree at most `2k`, and the
fresh-singleton extension are checked. Every target between `k+1` and the largest
computed value is attained in each tested prefix. For `k=2,3,4`, this checks
respectively `t=1,...,6`, `t=1,...,5`, and `t=1,...,4`.

These are exact finite checks, not a proof for unbounded `k` and `t`.
The prospective uniform argument uses bounded neighborhood capacity and a least
crossing of an integer sequence with upward increments at most one. It does not
assume monotonicity. No Lean proof is part of this experiment.
