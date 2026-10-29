# Cambie–Wanless bounded-array avoidance

Exact finite checks supporting a preregistration of Conjecture 5 in Stijn Cambie and Ian M. Wanless, *Avoiding a partial Latin square with an isotope of a given Latin square*, [arXiv:2610.20008v1](https://arxiv.org/abs/2610.20008v1), Sections 1–4.

Run with Python 3, without dependencies:

```sh
python3 check.py
```

The stored output is `check.out`. All arithmetic is integer or rational (`fractions.Fraction`); the program uses no floating point or randomness.

The fidelity checks reproduce the source's displayed order-3 partial Latin square, which none of the 12 Latin squares avoids, and its assertion that Conjecture 4 holds through order 3 (148 ordered pairs of squares at orders 2 and 3). These are published results, separate from the proposed Conjecture-5 resolution.

The checks supporting the expected proof enumerate all bipartite graphs through order 3 for the rook inclusion-exclusion identity, Bonferroni bounds and ordered-tuple collision estimate. They enumerate all eligible three-graph families at order 2 for the bounded-overlap identity, and all bipartite graphs through order 4 for the sufficient degree/common-neighborhood Hall criterion. Structural checks use every Latin square through order 3, every set-valued array through order 2, and every partial Latin square at order 3. They test conflict-graph incidence, the matching-to-symbol-permutation equivalence, forbidden-graph degree and edge-multiplicity bounds, invariance of the three array bounds under column/symbol exchange, and the exact correspondence between principal isotopes and conjugated row/symbol isotopes.

These checks test the finite identities in the proposed proof. They do not prove the uniform asymptotic avoidance estimate, permutation concentration, or the existence of a single threshold N(t) valid for all admissible L and A. Small-order unavoidable arrays are permitted by the conjecture's “sufficiently large” quantifier. No Lean proof of the conjecture or its negation is included.
