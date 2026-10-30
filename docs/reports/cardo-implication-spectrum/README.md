# Native implication spectrum: exact finite checks

Source: Carles Cardó, *Probabilistic equational spectrum, primality and
approximation in finite algebras*, arXiv:2604.08446v1, Definitions 2.1 and 3.1
and the unnumbered question in Section 7's final paragraph.

The proposed complete answer is the set of dyadic rationals in [1/4, 1] for
the two-element algebra whose only basic operation is material implication.
There are no constant symbols. These experiments are finite evidence for a
planned general proof, not a proof of the assertion at unbounded arity.

Run with Python 3's standard library:

    python3 check.py

`output.txt` contains the resulting stdout. Every probability uses exact
integer counts and `fractions.Fraction`; there is no floating-point evidence.

The fidelity check reproduces the source's semilattice formula
1 - (2^p + 2^q - 2)/2^(p+q+r) in all 27 cases p,q in {1,2,3} and r in {0,1,2}.
The implication check computes the least closure of projections under the
native operation for one through four variables, so it includes term
functions of arbitrary term size in these arities. The closure has respectively
2, 6, 38 and 942 functions. It agrees with the family of truth tables having a
coordinate whose value 1 forces the output 1. All pairs have precisely the
agreement counts predicted by the proposed spectrum at each tested arity.

A separate syntax construction uses only variables and implication to
realize every guarded truth table through four variables. It also constructs
equations for every permitted agreement numerator through six variables:
98 equation pairs in total. Every constructed equation uses all its counted
variables; adding a missing variable x as (x -> x) -> t preserves t. A scalar
Boolean evaluator independently checks every constructed syntax node against
its integer truth table at every valuation. The endpoint equations are
x = (x -> y), with probability 1/4, and x = x, with probability 1.

The candidate is expected to require structural induction, guarded Boolean
synthesis, and a finite-set partition argument valid at arbitrary arity.
The finite checks are not proposed as positive computational Lean content.
