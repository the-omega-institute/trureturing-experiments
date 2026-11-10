# Lu–Li Question 8.6: degree-supported eigenvector necessity

`check.py` independently checks the final yes/no clause of Question 8.6 in Zhipeng Lu and Pengxiang Li, *Adjacency-degree algebras and spectral determination of graphs*, [arXiv:2607.21494v2](https://arxiv.org/abs/2607.21494v2). It reproduces Proposition 8.1 for both nine-vertex source seeds and checks their first two private-leaf coronas. The opening classification question is outside scope.

Run with Python 3.12 or newer:

```sh
python3 check.py
```

The program uses only the Python standard library. Linear algebra uses integers and `fractions.Fraction`, with no numerical rank or approximate eigenvalues. It computes the smallest jointly adjacency/degree invariant space containing the all-ones vector, and independently verifies invariance of the resulting basis. Color refinement runs to an unchanged partition and verifies equitability. For each seed it checks the stated stable cells, the six-dimensional walk space, and the source eigenvector.

For both coronas it checks exact equality with the two lifted copies of the previous cyclic space and stable partition, and the degree-polynomial leaf projection. It verifies pendant/twin reduction certificates for distance heredity. On each nine- and eighteen-vertex graph it additionally exhausts all nonempty vertex subsets, checking distances in every connected induced subgraph directly by breadth-first search.

For each corona degree class S, the output supplies row and column indices of an identity minor in the adjacency block from S to its complement. Every vector supported on S and satisfying Af = λf must lie in that block's kernel. The identity minor makes this kernel zero, for every real λ simultaneously; this is stronger than exclusion inside W and the orthogonal complement of the all-ones vector. The hidden two-dimensional plane of the first G1 corona is checked directly.

`check.out` is the exact program output, ending with `ALL_CHECKS_OK`. These checks provide finite evidence and certificates. They are not Lean proofs and do not by themselves prove the infinite corona lift. The proposed general argument is part of the preregistration draft.
