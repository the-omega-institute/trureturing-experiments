# Clique roots of connected K4-free graphs

Exact checks for Hossein Teimoori Faal, *Sums of Polynomials and Clique Roots*, arXiv:2112.09744v1, Conjecture 4.1. The proposed counterfamily consists of a complete tripartite graph with three parts of size `2m`, together with one new vertex adjacent to exactly `m` vertices of one part.

Run `python3 check.py` with Python 3 and SymPy 1.14.0. The expected output is stored in `check.out`; success ends with `ALL_CHECKS_OK`. All evidence uses integers, rational numbers, symbolic algebra, or exact Sturm root counts.

The fidelity checks reproduce Proposition 3.3 on all connected labelled triangle-free graphs with at most six vertices, and the roots in Proposition 3.1 on paths with two through ten vertices. The family checks enumerate cliques directly from adjacency for `m = 2..8`, independently reject every four-vertex clique, count real roots exactly, and exhibit a disconnecting cut of size `m`. For `m = 2..5`, every deletion of fewer than `m` vertices is checked for connectivity.

Symbolic checks establish the coefficient, affine-substitution, derivative, critical-value, and discriminant identities in the parameter `m`. Finite graph checks do not prove the uniform graph-counting or connectivity statements. No Lean proof or exhaustive literature-clearance claim is supplied by this experiment.
