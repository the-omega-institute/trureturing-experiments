# Positive complete-label relation certificates

Run from this experiment repository with Python 3.12 or newer:

```sh
python3 rro-positive-certificates/exact_checks.py --output rro-positive-certificates/exact-results.json
```

The standard-library program loads existing literal-source helpers from the pinned
Git object `865f38257ead1b4e857ff915c23b0b04276bd7f3:rro-structural-recovery/exact_checks.py`
in memory. It leaves sparse checkout rules and the older experiment files intact.
The result records that supplier's SHA256. No third-party dependency is required.

The public mathematical source is
`docs/develop/theory/RECURSIVE_RELATIONAL_OBSERVATION_POSITIVE_RELATION_CERTIFICATES.md`
in trureturing. It proves uniform soundness over full injections. This program
checks the literal B0/B1 graphs, all 17496 endpoint incidences and their complete
capacity lifts, 576/324 positive witnesses and reconstructed 1944-edge relations,
all distinct left triples for the K3,10 boundary, ambient-slack counterexamples,
weakened capacity evidence, full-support and double-star failures, integer
occupancy counts, and actual original-root Read scans with retained indexed
partial/completed/unread domains. The label-only completion routines are sufficient
witness-pattern tests on supported positive histories, not a fast implementation
of the theory's exhaustive recognition test on every possible certificate.

The exact occupancy results refer to one selected witness's coverage event and
supply lower bounds for acceptance; they are not exact acceptance probabilities.
The symbolic proof, not random label tests, establishes full relabeling soundness.
No minimum-size, efficient-discovery, new Lean certification or originality claim
is attached to these finite checks.
