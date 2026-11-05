# Tree-beta to tree-MY corner scaling: exact checks

Source: Gyotoku, Sasada and Wesołowski, *Independence properties for tree beta
models*, arXiv:2609.15199v1, especially Lemmas 2.3–2.5 and the open question at
the end of Section 1.3. Source URL: https://arxiv.org/abs/2609.15199v1.

Run `python3 check.py` with Python 3 and SymPy. The stored result is `check.out`.
The program uses exact integer, polynomial and rational-function arithmetic;
it makes no floating-point comparisons.

The checks cover every labeled tree on 2, 3 and 4 vertices and every root:
20 trees, 75 rooted cases, and 284 induced principal subsets. Source fidelity
checks verify the inverse-map recursion and Delta factorization. Scaling checks
verify the leading coefficient of every induced Delta polynomial under
q_e = 1 − epsilon² c_e² and k_v = 1 − epsilon y_v. That coefficient is the
corresponding principal determinant. They also verify the limiting Schur
transport, its unit Jacobian, its rate identity, the exact diagonal Jacobian
before taking the limit, and two symbolic logarithmic kernel limits.

All checks pass. These finite identities are sanity checks for a proposed
general proof. They do not establish the universal analytic limit, dominated
convergence, integrability, or the absence of a prior settlement. The planned
positive settlement uses a general argument, not finite enumeration.
