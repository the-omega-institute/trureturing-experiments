# Atkins: minimum algebraic representation of 30_65

Independent exact checks supporting a preregistration of Eli Atkins,
*Relation Algebra Representations from Distance-Regular Graphs*,
[arXiv:2605.14190v1](https://arxiv.org/abs/2605.14190v1), Question 7.5.

Run with Python 3, using only its standard library:

```sh
python3 check.py
```

The expected output is in [output.txt](output.txt); success ends with
`ALL_CHECKS_OK`. Every evidential calculation uses integers, sets, or exact
rational arithmetic through `fractions.Fraction`.

The program constructs the Hoffman–Singleton graph from five pentagons and
five pentagrams, with cross edges defined modulo five. It independently
checks its strongly regular parameters (50,7,0,1), removes a vertex and its
seven neighbours, and checks the resulting 42-vertex graph's intersection
array {6,5,1;1,1,6}. Distance colours reproduce all 28,224 ordered-base/product
conditions of the full 30_65 composition table, including identity products
and the required witnesses. It verifies the seven six-point c-cliques,
the b-matchings, and the adjacency quadratic on the fibre-zero-sum subspace.
Three explicitly computed and fully verified permutations generate the four
ordered-pair colour orbits, so the source's algebraicity assertion is also
checked directly.

For the proposed lower bound, the program checks the entire small-clique
parameter list obtained from the general spectral argument: r=3 has no
permitted t; r=4 and r=5 fail either integral common-neighbour counts or
integral eigenvalue multiplicities. A separate finite parameter loop checks
all count-feasible possibilities below 42. This loop does not prove that
arbitrary algebraic representations satisfy those restrictions. That
structural implication, the spectral argument, and their unbounded scope
belong to the proposed mathematical proof. The attained parameters are
r=k=6, mu=1, with restricted multiplicities 21 and 14.

The spectral restrictions are established prior work: Coutinho, Godsil,
Shirazi and Zhan, *Equiangular Lines and Covers of the Complete Graph*,
[arXiv:1504.00085](https://arxiv.org/abs/1504.00085), Theorem 3.1;
[DOI:10.1016/j.laa.2015.09.029](https://doi.org/10.1016/j.laa.2015.09.029).
The proposed new deduction connects Atkins's literal orbital definition and
composition table to the universal 42-point minimum. No Lean proof is
contained in this report.

Producer: independent codex-cli research seat, flight
546493e74ee936eea59eca96; prior disclosure repo-prior-exposed. The program is
newly written here and embeds no copied third-party code or graph data.
