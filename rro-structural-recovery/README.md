# Exact evidence for the RRO source-capacity theorem

The mathematical source is `RECURSIVE_RELATIONAL_OBSERVATION_UNIFIED_JOINT_LAW_BRIDGE.md`, §§15–18, with the literal source definitions UJB.7502 and the shared full-injection prior UJB.9001. Coordinates in the programs are zero-based and labels are binary integers with the most significant coordinate first.

`exact_evaluator.py` implements the unrestricted finite indexed-cache formula: exact sample equality partitions, source geometric type assignments with remaining capacities, and the subset recurrence for compatible distinct labels. It accepts a JSON cache and an optional rational prior. Its loops may be enormous; it makes no practical speed or storage claim. It never enumerates full hidden injections.

`exact_checks.py` independently constructs both literal graphs from their bilinear formulas. It checks geometric keys against equal neighborhoods, left/right common-neighbor counts, capacities, positive and nonforcing witnesses, and the mixed predictions. For small graphs it compares the capacity expansion against graph-homomorphism partition inversion. For labels it compares the subset recurrence against collision-partition inclusion-exclusion on the actual 128/256 alphabets. Its at-most-four-sample inversion evaluator is a finite verification route, distinct from the unrestricted evaluator and from the ordinary proof for all masks. Complete and partial four-sample normalizers are compared using both routes. Independent uniform-image bucket sums cover the joint-page obstruction; literal complete maps cover the old-slot obstruction. Recursive term readout is compared with the full-menu operand classifier, including absent and noncanonical suffixes.

Run with Python 3.12 or newer, using only its standard library:

```sh
PYTHONDONTWRITEBYTECODE=1 python3 exact_checks.py --output exact-results.json
```

`exact-results.json` contains the exact useful counts, rational outputs, finite coverage and environment. These finite checks support, but do not replace, the ordinary proofs for arbitrary finite caches, all prior parameters and the full continued-use interface. No formal compilation, polynomial bound, runtime improvement, optimal total memory, world novelty, native randomness generation, hidden-realization inversion or physical-spacetime inference follows from these data.
