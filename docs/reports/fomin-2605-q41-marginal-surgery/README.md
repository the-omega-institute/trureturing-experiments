# Fomin Question 4.1: marginal surgery

Source: Dmitri V. Fomin, *Closed polylines with fixed self-intersection index*,
arXiv:2605.05506v3, Question 4.1 and construction K9.
The arXiv version history identifies v3 as revised 17 May 2026.

Run with Python 3 and its standard library:

```sh
python3 docs/reports/fomin-2605-q41-marginal-surgery/check.py
```

The program writes `result.json` and ends with `ALL_CHECKS_OK`.
`check.out` contains its complete stdout. All calculations use integers and
`fractions.Fraction`. Two independent intersection tests, strict orientation
products and exact segment parameters, must agree for every edge pair.

The fidelity check reads the integer vertex order from the source's K10 figure
and reproduces its stated type (12, 4): 24 proper crossings, four on every edge.
The refutation check uses ten explicitly listed integer vertices. It checks one
closed cycle, distinct vertices, no three collinear vertices, no nonincident
vertex on an edge, no parallel edge directions, 25 distinct proper crossings,
and exactly five crossings on every edge. Its complete crossing orders and
exact parameters are included in the JSON.

Exactly five crossings satisfy K9/1: edge pairs (0, 5), (1, 6), (2, 7),
(3, 8), (4, 9). Each is first on both edges. For every candidate, the program
removes AB and MN and inserts AN and BM using the endpoints B and N adjacent to
the crossing. Each replacement leaves two cycles with five vertices each.
The same surgery accepts an opposite-orientation positive control. The K10
fidelity example also has seven crossings passing both conditions, so this
implementation does not reject connectivity universally.

This is an exact finite certificate for a proposed refutation. It is not a Lean
proof, a claim of literature novelty, a minimality assertion, or verification
of the proposed infinite doubled-star family. The family belongs to the
preregistration's mathematical mechanism discussion.
