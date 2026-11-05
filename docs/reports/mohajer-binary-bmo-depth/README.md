# Binary-tree BMO: first amplification depth

Source: Mohammadreza Mohajer, *Sharp Constants and Extremizers in Homogeneous Tree BMO*, arXiv:2610.14459v1, Section 6(2), with leaf formulas in Section 2.3 and the depth-four witness in Section 4.4.

Run `python3 check.py` with Python 3. The standard library is sufficient. The stored output is `check.out`; success ends with `ALL_CHECKS_OK`. All numerical evidence uses integers or `fractions.Fraction`.

The program reproduces every input-node oscillation of the published depth-four vector and its exact amplification ratio 41/40. For the proposed universal depth-three contraction argument, it expands all six nonnegative linear certificates, checks coefficient domination on every complementary sign cone, verifies the 128 tree automorphisms and the six stabilizer orbits, covers all 168 marked triples, checks the opposite-half two-tail certificates, and verifies the sharpness vectors. It also checks source tail identities and contraction on all 6,561 ternary eight-leaf assignments and 4,096 deterministic rational assignments, including ties and mixed tails of sizes two and three.

The coefficient and orbit checks support a general argument over arbitrary real leaf values. The sample checks alone do not establish that universal statement. No Lean proof or kernel verification of the proposed settlement is supplied here. The exact depth-four constant and all-depth constant remain outside this experiment.
