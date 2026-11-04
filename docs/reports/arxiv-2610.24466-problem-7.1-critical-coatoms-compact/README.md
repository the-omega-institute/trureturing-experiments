# Critical coalition coatoms: exact finite checks

Source: Shanxia Wang and Yuqiang Dou, *The Lattice of Coalition Logics between MCL and CL: Coatoms and Compact Meets*, arXiv:2610.24466v1, Problem 7.1 and Table 1. Primary text: https://arxiv.org/html/2610.24466v1 .

Run `python3 check.py` with Python 3.9 or newer; only the standard library is used. The stored output is `check.out`. Every calculation uses integers, Boolean bitsets, and exhaustive finite sets.

The independent full-profile implementation reproduces Table 1 and the Ser/IA/Det validity types for all admissible combinations of primes 2, 3, 5, 7 and agent counts 2, 3, 4. It includes every coalition and every dummy-agent action. These 27 models give 14,000 operation-table comparisons and 7,593,612 classical-axiom valuations.

For the proposed affirmative settlement, the program checks the cyclic-shift finite basis and every centrality/guard equation, including disjoint guard variable tuples, on all these models and 54 direct products. It also exhausts 66,049 possible tuples satisfying necessary MCL conditions on the four-element Boolean algebra with two agents. Exactly 784 satisfy MCL. For prime 2 and the IA family, exactly 122 satisfy the proposed equations, agreeing exactly with the union of 121 CL models and two cyclic-basis models (one model belongs to both).

These checks support the proposed finite-basis construction. They prove neither compactness of the logic lattice nor the uniform result over all primes and finite agent sets. The general proof requires Boolean normal forms, central-factor decompositions, subdirect representation, and the structural-logic/variety correspondence. No Lean proof of the claim or its negation is included.
