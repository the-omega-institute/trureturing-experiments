# Exact local GW/fGW certificates

The check supports the exponent-extension proposal for Section 6, first bullet,
of Hohmeier, Fraiman and Moosmüller, *k-Nearest Neighbors in Gromov–Wasserstein
Space*, arXiv:2606.10295v1. It uses only Python's standard library and exact
integer and rational arithmetic.

Run `python3 docs/reports/hfm2606knn/check.py` from the repository root. The stored
output is `check.out`; success ends with `ALL_CHECKS_OK`.

For n from 1 through 4, the program enumerates every labelled metric with
independent off-diagonal entries in {1,2}, and compares it with its uniform
positive rational perturbation. All n! permutation couplings and every ordered
pair of these permutations are evaluated. GW exponents are 1,2,3,4. Fused pairs
are (1,1),(2,1),(2,2),(4,2) with alpha=1/2 and one-dimensional rational features,
and (1,2) with alpha=1/4 and zero features. In the last case square roots are
rational and checked by integer square roots, without numerical approximation.

Each example obtains a global certificate for the minimum over the continuous
coupling polytope. If the best diagonal cross-cost is m and all off-diagonal
cross-costs exceed m, every real Birkhoff convex combination has cost at least m;
a minimizing permutation attains m. The program verifies the stronger uniform
cross-cost gap from the local-equality argument. This is an exact certified
solver for these local instances; it is not a solver for arbitrary GW inputs.
It also independently compares direct four-index coupling sums with their
Birkhoff cross-term expansions at rational mixtures.

Fidelity checks reproduce the paper's GW_2/permutation comparison and its
nonnegative-cross-term factorial bound on these instances. An additional complete
two-point calculation checks all edge lengths a,b in {1,...,5}, obtaining
GW_2 squared = (a-b)^2/2, with the minimum zero precisely for isometric supports.

The check verifies 675 local instances, 320148 ordered cross-cost certificates,
2673 direct objective expansions, and 100 source fidelity checks. It does not
verify arbitrary real exponents, differentiation, or a classifier risk limit.
Those require the general mathematical argument and subsequent Lean work.
