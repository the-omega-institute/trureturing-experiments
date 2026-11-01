# Yang–Brandenburger agreement without certainty reflection

Run with Python 3, using only the standard library:

```sh
python3 docs/reports/yang-brandenburger-agreement-without-reflection/check.py
```

The program uses exact integers and `fractions.Fraction`. The saved output is
`check.out`, whose final line is `ALL_CHECKS_OK`.

The fidelity check reproduces Example 1 in Section 2 of Erya Yang and Adam
Brandenburger, *Conditional Probability Spaces and the Structure of Agreement*,
arXiv:2605.30017v2: the terminal sets are {b,c} and {a,b,c}, giving common
certainty of probabilities 1 and 0 when 1-closedness is absent.

The candidate check enumerates every preorder on 1, 2 and 3 labelled states,
every surjective rank assignment into an initial finite interval, and the two
positive weight vectors (1,…,1) and (1,…,n). On each conditioning event, the
lowest-rank states receive probabilities proportional to their positive weights.
It retains 1-closed CPSs and deduplicates their probability tables. It validates
all CPS axioms, the source's truth-plus-certainty and overlapping-atom averaging
lemmas, the proposed common-support and local-saturation steps, and agreement
under local consistency. There are 1, 12 and 151 distinct retained models, with
218,122 event/state tests and no failures. Probability values at the designated
state suffice: any common-certainty assessment must already equal those values.
This sample does not exhaust arbitrary real CPSs and is not a universal proof.

A separate four-state uniform example satisfies the retained hypotheses with
common posterior 1/2 and fails certainty reflection for both agents.

The program also isolates a source discrepancy: the displayed Bob family in
Example 2, Section 3, is not 1-closed. Its measure on Ω gives {a,c,d} probability
one, but that event is missing from the conditioning family. This observation
concerns that example's asserted hypothesis; it is not a counterexample to the
reflection-removal question, which retains 1-closedness.

The intended settlement of Section 5.3's first question is a general argument
for arbitrary finite state spaces. No Lean proof or probe is included here.
