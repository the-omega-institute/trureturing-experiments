# Two exact hexagonal-growth checks

Run `python3 docs/reports/a400701-hexgrowth/check.py` with Python 3 (standard library only). The stored `check.out` is its complete output. Success requires exit 0 and the final line `ALL_CHECKS_OK`.

The native simultaneous threshold-two update reproduces all 51 printed terms of [OEIS A400701, revision 30](https://oeis.org/A400701) and all 52 printed terms of [OEIS A374338, revision 15](https://oeis.org/A374338). A separate scan checks 24 native updates. Full occupied sets match the proposed six-inequality shapes at both seeds for every horizon 2 through 100. The cross-sequence identity holds at all 99 horizons 2 through 100; all three A374338 residue formulas hold through horizon 100.

The circle check computes the complete intersections of equal unit circles at all 18 nonzero triangular-lattice displacements of squared length at most 4, using `fractions.Fraction`. It then compares eight full geometric circle iterations with the native two-site update. A separate congruence-class enumeration checks the square-minus-corners cardinality formula through horizon 100. All arithmetic used as evidence is exact; no floating point is used.

These are finite computations, not a proof of either unbounded conjecture. The Euclidean-circle interpretation retains old vertices and all previously drawn circles. The source's iteration-2 illustration explicitly labels the count as eight total vertices. Only OEIS numerical data and the rule are used; neither proposed generating function is assumed.
