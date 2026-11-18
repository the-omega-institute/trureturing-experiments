# Exact checks for Tărnăuceanu Conjecture 4

Source: Marius Tărnăuceanu, *A note on subgroup commutativity degrees of finite groups*, arXiv:1802.03552v1, Conjecture 4 and the preceding Example, equation (2).

Run `python3 check.py`. Python 3.12 or later; standard library only. The output is in `output.txt`. All arithmetic used as evidence is integer or rational arithmetic.

The program enumerates subgroup lattices, set products, normal subgroups and quotient groups from multiplication tables. It independently reproduces the source's `sd(D8) = 23/25` and `sd(A4) = 16/25`. It checks the scalar-group formula and section minima for the groups of orders 12, 56 and 80, their centers, and their proper subgroups. It checks coprime subgroup splitting and multiplicativity on A4 × C5.

For the proposed large witness, the program verifies the Gaussian subgroup counts using both a recurrence and the product formula, multiplicative orders, disjoint prime supports, and the exact rational inequality above 23/25. These checks use the mathematically classified subgroup lattice; they do not enumerate the large groups. Finite samples of the proposed unbounded family certify the graph-subspace lower bound for several rational thresholds and numbers of factors. A finite program does not prove the universal existence theorem, the full group classification, or Conjecture 4's negation. Those remain obligations of the subsequent Lean probe and implementation.

The intended construction is a coprime product of centerless irreducible scalar Frobenius groups, each a Schmidt group. The general argument, rather than finite positive enumeration, supplies the proposed refutation family.
