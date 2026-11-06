# Exact finite checks for the proposed Kalman extension

Source: J. Wang, J. Shi and M. Wang, arXiv:2609.09529v1, Definitions 2.1, 2.3, 3.8, 4.1 and 4.8, Theorem 4.17, Open Problem 6.1.

Run with Python 3: `python3 check.py`. The standard library suffices. Every operation and comparison uses integers or finite tuples of integers; there is no floating-point evidence. The program and its output are independent checks, not a proof or a formal settlement.

The enumerated families are all commutative integral residuated monoid operations on chains of cardinality at most five, all FARL unary expansions of those reducts, all FARL unary expansions of the four-element Boolean lattice with meet multiplication, and all FARL unary expansions of the following six-element nondistributive reduct. Its lattice is 0 < p,q,r < e < 1, where p,q,r are incomparable; multiplication of two nonunits is 0 and 1 is the unit. Residuation is computed by its defining adjunction. Enumeration of the unary maps is exhaustive within these families: g(0)=0, x≤g(x), and preservation of joins are necessary for every source expansion, and the unique right adjoint f is then computed and all five FARL axioms checked.

For every accepted expansion the program checks all source axioms, the twist operations, all centered target axioms, Leibniz, D, S and CK, the positive-cone source axioms, and preservation of every operation and bijectivity of rho and sigma. It also checks the multiplication bridge used in recovering sigma. There are 124 expansions across the enumerated families; the three nondistributive expansions each have a 27-element twist.

For source fidelity the program reproduces the displayed Łukasiewicz-square example of Section 4.1: coordinates 0,1,2 represent 0,1/2,1, f is the diagonal coordinate minimum, and its left adjoint is the diagonal maximum. It directly checks the displayed residual formula, double negation, all three FAIRL equations, the derived g, all five FARL equations and the two recovery isomorphisms, with a 36-element twist. This is an additional test, outside the 124-expansion enumeration.

The finite tests do not enumerate all finite lattices, all target-category objects, or all homomorphisms, and do not establish the equivalence for infinite algebras. The proposed positive settlement requires a general mathematical argument. No third-party program source was read, copied or incorporated in this check.
