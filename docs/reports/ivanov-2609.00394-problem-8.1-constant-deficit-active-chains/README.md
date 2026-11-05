# Constant-deficit active pseudosimilarity chains

Exact finite checks for Sergey Ivanov, *Graphs with Long Pseudosimilarity Chains under Consecutive Vertex Deletions*, [arXiv:2609.00394v1](https://arxiv.org/abs/2609.00394v1), Problem 8.1. Source version submitted 31 August 2026, licensed CC BY 4.0. The source's Section 3, Example 3.1 supplies the fidelity target.

Run with Python 3 and NetworkX 3.6.1:

```sh
python3 -m venv /tmp/ivanov-check-env
/tmp/ivanov-check-env/bin/pip install networkx==3.6.1
/tmp/ivanov-check-env/bin/python check.py > check.out
```

The program uses only integer labels, exact degrees and exact adjacency comparisons. NetworkX's VF2 enumerator independently checks automorphism groups and deletion-card isomorphisms. Its output ends with `ALL_CHECKS_OK` after all assertions succeed; an assertion failure produces a nonzero exit.

The six-clock Q has vertices a_i = i and b_i = 6+i, with indices modulo 6. A induces a six-cycle, B is independent, and b_i has neighbors a_i, a_(i+1), a_(i+3). K is its complement. The proposed family adjoins a path x_0,...,x_(m-1), represented by integers 12,...,12+m-1, and attaches x_i to a_(i mod 6).

Checks and results in `check.out`:

- Fidelity: the source's 17-vertex, 29-edge example, its four specified active deletions and the asymmetric connected stages.
- Core: all 720 permutations of A, with the B action forced by distinct neighborhoods; exactly the six simultaneous rotations preserve K. VF2 independently enumerates the same six automorphisms.
- Family: m = 1,...,100 (orders 13,...,112), connectivity, degree separation, trivial automorphism groups and, when m ≥ 2, the explicit endpoint-card isomorphism and an independent VF2 isomorphism test.
- Entire chains for starting orders 14,...,18: every intermediate graph, including the terminal graph, is connected and asymmetric; every deleted endpoint has the surviving x_0 as a pseudosimilar mate. Chain length is n-13.

These finite results corroborate the proposed construction. They do not prove a statement over unbounded graph orders. The intended mathematical argument uses a symbolic automorphism classification and a uniform path-deletion induction.
