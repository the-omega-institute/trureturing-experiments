# Elphick's graph-energy Conjecture 6.2

Source: Felix Goldberg, *Graph energy estimates via the Chebyshev functional*,
arXiv:1407.7430v3, Section 6, Conjecture 6.2, https://arxiv.org/abs/1407.7430v3.
The conjecture is E <= 2m/sqrt(lambda_1) for connected ordinary simple graphs.
Energy is the sum of absolute adjacency eigenvalues, with multiplicities.

Run `python3 check.py` with Python 3 and SymPy (output generated with SymPy 1.14.0).
All evidence uses integers or exact symbolic algebra. `output.txt` is the output.
The final line is `ALL_CHECKS_OK`.

The fidelity check reproduces the complete-bipartite energy and equality in the
source's lower bound, on 20 graphs, as well as trace(A^2)=2m. The source's adjacent
claim `m=2pq` is a typo: the actual unordered edge count of K_(p,q) is pq. The
check uses the actual adjacency matrix and pq edges.

For T_t, subdivide every edge of the t-arm star once. Exact characteristic
polynomials are checked at t=1..12,15,16,20. At t=15 an independent rational
31-vector eigenbasis verifies rho=4, E=36, m=30, and proposed upper bound 30.
The symbolic quotient and gap identities support the proposed unbounded family
t=k^4-1: E*sqrt(rho)/(2m) = k/2 + k/(2(k^2+1)).

Adding a leaf at the center gives U_t, a nonsingular family. Exact characteristic
polynomials and nonzero determinants are checked at t=1..8,16,20; exact radical
bounds are checked at t=16,17,20,100,10000. These are experiment checks and
algebraic identities, not Lean proofs of a universal family assertion. They do
not reproduce the source's exhaustive eight-vertex graph computation.

The source is cited, not redistributed. The check program is independently
written for this report; its use follows the experiments repository license.
