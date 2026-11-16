# Disruption-recognition compiler checks

Source: Jose Pascual Gumbau Mezquita, *Rice's Theorem under Self-Modification:
Elevation Operators and a Normal Form*, arXiv:2609.11326v3, Open Problem 9.1.

Run `python3 docs/reports/gm2609op91/check.py` from the repository root.
The program uses only Python's standard library and exact integer operations.
`check.out` contains its output.

The fidelity check reproduces the source's identity-property example: an inert
wrapper computes identity; its instrumented rewrite computes identity precisely
when the embedded event halts, and otherwise is nowhere defined.

The independent checks use seven explicitly specified register-machine events,
including halting events and exact cycles. They test every nontrivial subset of
five distinct computable behaviours and every accepted base in each subset.
For each property the instrumented rewrite has a fixed target behaviour selected
according to whether the nowhere-defined behaviour belongs to the property.

Compiler A returns a program that waits for a selected event and then runs the
instrumented output. Compiler B returns a constant-zero program when the event
halts within the numerical input bound, and otherwise the instrumented output.
The checks compare both compiled programs against independently calculated
outputs and check Compiler B's finite-exception/constant-tail form.

These finite operational checks support the proposed compiler formulas. They
prove no undecidability statement and decide no arbitrary behavioural property.
In the general argument, ruling out every wrapper uses the non-recursiveness
of the halting set and the fact that every finite modification of a constant
predicate is recursive; it does not follow from enumerating the fixtures.
