# Exact checks for Fish's finite-ring conjecture

Source: Alexander Fish, *On product of difference sets for sets of positive density*, arXiv:1702.02544v2, §4, Conjecture 1. The conjecture asks for a subring of uniformly bounded index inside the set of individual products of two differences of a dense set.

Run `python3 docs/reports/fish1702c1/check.py` from the repository root. The program uses only integer arithmetic and `fractions.Fraction` from Python's standard library. The output is `check.out`; its last line is `ALL_CHECKS_OK`.

The source-fidelity check reproduces §3's identity between positive-measure cyclic return times and the periodic lift of the difference set in 33,798 cases. It also checks the high-density specialization of Corollary 1.2 on 462 cyclic subsets.

The proposed refutation uses Boolean rings with xor addition and support-intersection multiplication. The program checks every Hamming-ball radius in dimensions zero through eight, directly checks the four-variable product in dimensions at most five, enumerates every additive subspace and every Boolean subring in dimensions at most four, and verifies the exact maximal subring cardinality for each support bound. Exact binomial calculations check the stated family for parameters two through eight; recurrence and induction-polynomial checks cover parameters zero through 512.

The separate density-one-half family uses the multiplication

`(a,v,c)(b,w,h) = (ab, aw+bv, ah+bc+v·w)`

over the two-element field. For dimensions one through four, the program exhaustively checks associativity, distributivity, commutativity, the identity, units, the difference product, radical triple products, and the exact maximal contained-subring cardinality.

These are finite checks supporting the preregistered ordinary argument. They do not prove the conjecture's negation in Lean, establish the infinite-parameter result by enumeration, or certify the completeness of a literature search. The source-specific multiplicative reduction and exact subring-size argument are distinguished from the classical additive Hamming-ball construction.
