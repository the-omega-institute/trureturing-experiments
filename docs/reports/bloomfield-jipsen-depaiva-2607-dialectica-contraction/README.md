# Dialectica contraction: exact checks

Source: Colin Bloomfield, Peter Jipsen and Valeria de Paiva, *Dialectica Categories over Heyting Algebras*, arXiv:2607.23430v1, Section 4.

Run `python3 check.py`. The standard-library program uses finite sets and integers throughout. The recorded output is `check.out`; success requires exit 0 and the final line `ALL_CHECKS_OK`.

The fidelity checks reproduce the displayed tensor table of D′(2), the ten-element nondistributive lattice D′(3), and the source's membership-dependent Set contraction selector for all 104 relations with witness cardinality at most 2 and counterobject cardinality at most 3.

For the proposed counterexample in the covariant diagram category Set^J, where J is the span r → s and r → t, the constant counterobject has underlying set {0,1,2}; the three relation components are {2}, {0,2}, and {1,2}. Identity transition maps force one shared selector. Exhaustion of all 19,683 possible shared functions gives no contraction selector. Each stage separately admits a selector; the independently counted stage totals are 768, 81, and 81.

A second calculation computes the permitted output intersection at every input pair. It rules out a selector for all 64 block-size triples (p,q,r) in {1,2,3,4}³. The obstruction at a mixed pair from the first two blocks is that the two leaves force disjoint output blocks. All three relation components are nonempty and proper.

These computations are evidence for preregistration, not a Lean proof. The arbitrary positive block-size family, the base category's universal properties, and the connection to Mathlib's subobject inequality remain obligations of the later proof. No proof of the source's assertion or its negation is included here.
