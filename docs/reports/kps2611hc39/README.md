# Honeycomb critical spacing 39: exact checks

Source: T. Kittipassorn, N. Phetmak and K. Sriamorn, *Covering the Square Lattice and the Honeycomb with Non-Overlapping Unit Disks*, arXiv:2611.16894v1 (12 November 2026), Section 7 Question (1), Lemma 3.14 and Appendix A Table 6 rows 69–70.

Run with Python 3.12 or newer, with no external dependencies:

```sh
python3 docs/reports/kps2611hc39/check.py
```

`check.out` is the exact output. Every numerical decision uses integers or `fractions.Fraction`; no floating-point result supplies evidence. The program was written independently from the mathematical definitions. The source's ancillary program and unlicensed third-party Lean code were not read or copied.

In coordinates a+bω, squared distance is Q(a,b)=a²−ab+b². Rescaling the question gives disk radius squared 39/4. A minimum enclosing disk of a nonempty finite point set has a support consisting of a singleton, a diameter pair, or three noncollinear points. Translate one support point by the honeycomb period lattice to +1 or −1. Every other point is at Q-distance at most 39 from that anchor. Since Q(a,b)≥3a²/4 and Q(a,b)≥3b²/4, both relative coordinates lie between −7 and 7. Enumerating these supports and counting their points gives a finite exhaustive reduction for the disk counts, including all disks of strictly smaller radius. This reduction remains a mathematical obligation for a later Lean implementation.

Results:

- The endpoint maximum is 26. Over every strictly smaller radius the maximum is 25; it is already attained at squared radius 9. The support census contains 1,310 circles over the two anchor classes.
- Every enumerated 26-point disk has an antipodal honeycomb pair. Every endpoint circumcircle triangle has a squared side length 39. The norm-26 alternative has no integer solutions, and Q(2z)=156 has 12 solutions, all with even coordinates.
- There are 12 norm-39 generators z and six eligible classes of doubled centers modulo twice the honeycomb period lattice. For all 72 combinations, the radius-squared-39/4 disk has 26 points and represents exactly 25 of the 26 honeycomb classes modulo zℤ[ω]. Explicit antipodal boundary points are constructed and checked in every case. The quotient is enumerated modulo 39, using the adjugate of the integer multiplication matrix of z; hence its 39 classes are checked without approximate geometry.
- The triangle-arithmetic step of the conditional prime extension is checked at p=7,13,19,31,37,43. No disk-count hypotheses are asserted at those other primes.

The checks support the proposed negative answer. They do not establish the infinite compactness/Voronoi equality reduction in Lemma 3.14, nor any Lean settlement. In particular, periodic center sets are not substituted for the arbitrary center sets quantified over in the original question: the reduction to a regular packing must be proved separately.
