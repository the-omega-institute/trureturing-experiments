# Countable closed spherical subsets and sequence l1

Source: Ian Doust and Anthony Weston, “Squaring the circle: embedding S¹ in ℓ₁”, arXiv:2608.10422v1, Section 5, question (2). The source was submitted on 11 August 2026. The problem asks about the great-circle metric, with target the real sequence space ℓ₁(N).

Run with Python 3 and the standard library:

```sh
python3 docs/reports/doust-weston-countable-spherical-l1/check.py
```

Every distance is represented by its exact rational coefficient of π. The program uses integers and `fractions.Fraction`; it uses no floating point or numerical quadrature. `output.txt` is the program's complete standard output.

- Fidelity: the explicit source embedding of the dyadic equatorial antipodal pairs together with the north pole preserves all 729 ordered distances for indices 0 through 12. The infinite geometric tails are summed exactly rather than discarded.
- Candidate identity: for every regular equatorial 2n-gon with both poles, for n from 1 through 16, the hemisphere-cut construction preserves all 7,136 ordered distances. Its normal-sphere cells are half-lunes of normalized area 1/(4n). The angular widths and all resulting distances are rational multiples of π. Three interior representatives per cell have identical sign patterns.
- Degenerate cases: all 64 subsets of the six-point octahedron, including empty, singleton and antipodal subsets, satisfy the identity after equal restricted sign patterns are combined; 672 ordered distances are checked.
- Negative control: halving the cut weights fails on an antipodal pair.

The proposed universal argument must justify a countable measurable partition by full sign patterns, its summable weights, and the angular separating-measure formula for arbitrary pairs. These finite tests do not prove that argument, measure zero of an infinite union, or the universal assertion. No Lean proof or kernel verification is contained here.
