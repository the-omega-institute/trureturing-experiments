# Quantum QBCE versus quantum individual sufficiency

`check.py` checks the proposed refutation of the general converse left unresolved
in Furkan Sezer, arXiv:2608.04973v1, Section 9.2. Run with Python 3 and SymPy 1.14:

```sh
python3 docs/reports/sezer-2026-qbce-qis-completeness/check.py
```

All calculations use exact integers, rationals, and symbolic expressions. The
stored output is `check.out`, ending in `ALL_CHECKS_OK`.

The fidelity checks reproduce the source's depolarizing-channel composition
(Section 10.2) and the trace identity in the forward direction of Theorem 6.2.
The counterexample check verifies PSD and trace-one ensembles, exact trace norms
6/5 and 8/5, symbolic two-ray cone identities, the pure-state block reduction,
and all player marginals for two, three, and four product factors. A signed
rational grid checks the coefficientwise positivity criterion in 245 cases.
The rational parametrization s = 2q/(1+q²), 0 < q < 1, covers the entire proposed
open parameter interval 0 < s < 1.

For n ≥ 2 and 0 < s < t < 1, use the two state-conditioned pure products of
|0⟩ and (√(1−s²), s), and the analogous target with parameter t. Pure-state
extremality restricts decision blocks to their ensemble rays. Local obedience
therefore imposes the classical correlated-equilibrium inequalities separately
in each payoff state, independently of s. A positive trace-preserving local map
would have to increase the trace norm of a Hermitian difference from 2s to 2t.

The program verifies the stated algebra and finite checks. The universal
all-games implication, positive-map contraction argument, and continuous-family
conclusion are a mathematical proof plan, awaiting Lean verification after
preregistration. No Lean proof or priority certification is claimed.
