# Liu–Diwan–Zhang forcing-number checks

Source: Qianqian Liu, Ajit A. Diwan and Heping Zhang, *Graphs with large maximum forcing number*, [arXiv:2512.22761v1](https://arxiv.org/abs/2512.22761v1), Problem 5.5.

Run with Python 3.10 or later, without third-party packages:

```sh
python3 docs/reports/liu-diwan-zhang-forcing-number/check.py
```

The program enumerates every Boolean bipartite adjacency matrix containing the identity matching for half-orders 1 through 4. Relabeling any perfect matching to the identity shows that these include a representative of every bipartite graph isomorphism type of those orders with a perfect matching. They are not an enumeration of every labeled matrix without this normalization.

Perfect matchings are enumerated as permutations. A forcing set is checked directly against every perfect matching, and every subset of the selected matching is considered in increasing cardinality. This reproduces Corollary 2.2 and Theorem 5.1 on all 4,165 graphs, with 16,774 matching forcing numbers. The script also reproduces Remark 5.6 for the four pairs (n,k)=(3,2),(4,3),(5,2),(7,3).

The proposed inequality `(f(G)+1)(n-F(G)) >= n` and every applicable Problem 5.5 implication pass. For one selected minimum forcing set in each graph, the program constructs the residual topological order and closed intervals for every perfect matching. All overlap counts, strict endpoint coloring decisions, induced uniqueness tests and complement-forcing tests pass: 16,774 matchings and 38,179 nonempty interval classes.

All evidence uses exact integers and finite combinatorics. The output is in [output.txt](output.txt), whose final line is `ALL_CHECKS_OK`. These are finite sanity checks for a proposed general argument; they do not prove the universal statement. No Lean proof is part of this report.
