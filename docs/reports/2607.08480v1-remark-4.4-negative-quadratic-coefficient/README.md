# Negative quadratic coefficient in the cubic independence quotient

Source: O. P. Bhardwaj, T. Chau, S. Ikram, G. Lather, V. Nanjangud and C. Venugopal, *Multiplicity of negative one of independence polynomials of graphs*, [arXiv:2607.08480v1](https://arxiv.org/abs/2607.08480v1), Remark 4.4.

Run `python3 check.py` with Python 3 and SymPy 1.14.0. The calculation uses integers and exact symbolic algebra. Its captured output is `check.out`; successful completion ends with `ALL_CHECKS_OK`.

The source fidelity checks reproduce the star independence polynomial for orders 2 through 10 and Theorem 4.1 on all 1,099 labelled simple graphs of orders 1 through 5 (870 satisfy its multiplicity hypothesis).

The candidate family is the join of a clique of size `6(t-1)^2` with the disjoint union of cliques of sizes `3,4,t,t`, for every integer `t >= 5`. The program checks the polynomial identity symbolically. Substitution `t=u+5` gives polynomials with strictly positive coefficients for the negatives of the quadratic quotient coefficient and of the quotient's value at minus one. These exact certificates support the stated infinite parameter range; the script does not infer it from sampling.

For `t=5`, explicit adjacency and an exhaustive independent-set recurrence give `(1,113,107,295,300)`. A separate enumeration of clique choices and the clique product/join formula give the same sequence. The graph has 113 vertices, 6,221 edges, independence number 4, and is connected without isolated vertices. Its quotient is `1+112x-5x^2+300x^3`, which evaluates to `-416` at minus one; the exact multiplicity is one. The script also checks the source's formulas for the cubic coefficients.

The strict coefficient valley supports a separate non-unimodality observation for the dossier's triage. That observation is independent of the positivity claim. Levit and Mandrescu's earlier joining-clique example, with joining size 25, is credited and checked: its value at minus one is 71, so it does not meet the required root constraint.

These are exact computational and algebraic checks for a preregistration draft. No Lean proof of the claim or its negation is included.
