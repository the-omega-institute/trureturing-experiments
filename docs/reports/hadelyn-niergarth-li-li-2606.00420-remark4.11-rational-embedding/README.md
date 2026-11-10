# Rational embedding obstruction for Remark 4.11

Source: Vixail Hadelyn, Harper Niergarth, Weiyou Li and Wenhui Li,
*Counterexamples regarding elementary symmetric partitions*,
[arXiv:2606.00420v3](https://arxiv.org/html/2606.00420v3),
Proposition 2.6, Definition 3.1, Proposition 3.4, Example 3.6 and Remark 4.11.

Run with Python 3.12 or newer:

```sh
python3 check.py
```

The program uses only integer arithmetic and `fractions.Fraction`.
`check.out` is its complete output. A successful run ends with `ALL_CHECKS_OK`.

The fidelity checks reproduce the normalization example with Z = 1680, the
edge calculation giving 7/33 and 77/75, and the equal-size triple-product
collision between (90,50,45,25,12,12) and (75,75,36,20,18,10).

The proposed general obstruction is that the two contributions from each
rational edge have a rational-square product. Common scaling preserves this
occurrence-level pairing. The target's square classes are (10,2,5,1,3,3),
so no pairing is possible. The program checks the edge identity on all 3,969
ordered pairs of reduced rational coordinates with denominators at most 10,
including degenerate equal coordinates and zero, and seven positive scales.
It enumerates all 4,166 loopless directed graphs on prefixes of a four-vertex
grid, checks all 166 Eulerian cases, and checks all 15 perfect pairings of the
target. It also checks 63 positive rational scales of the target and its first
200 positive integral multiples.

These finite checks support the preregistered general argument; they are not
a Lean proof of nonexistence or a classification of rational embeddings.
The source's two repeated parts 12 are separate occurrences. Antiparallel
edges, isolated vertices and edgeless graphs are permitted by its literal
directed-graph and Eulerian definitions. No connectedness condition is added.

The code is independently written from the paper's mathematical formulas.
It contains no third-party program source. SPDX-License-Identifier: Apache-2.0.
