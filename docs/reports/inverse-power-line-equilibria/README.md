# Exact checks for inverse-power line equilibria

The target is the unrestricted extension to every real exponent s > 2 of
Alper Ferudun, *Equilibria of Inverse-Square Repulsion on the Line Are
Arithmetic Progressions*, arXiv:2610.06944v1, Remark 5.3.
The source proves the classification for 1 < s ≤ 2. Its Lemma 2.3 states
the inverse-square conductance identity reproduced here.

Run with Python 3, without external dependencies:

```sh
python3 docs/reports/inverse-power-line-equilibria/check.py
```

`output.txt` contains the exact output. All arithmetic used as evidence is
integer or `fractions.Fraction` arithmetic. The program checks:

- the source conductance formula on 81 positive rational distance pairs,
  including equal distances;
- the aligned, truncated force-difference identity and its reciprocal-gap
  transform on 8,190 cases (all gap lists of length 2–5 with entries 1–3,
  every vertex, and integer exponents 2–6);
- exact cancellation for 3,780 symmetric arithmetic-progression windows;
- the induced ground-state energy identity on 6,003 cases (all graphs
  on at most three vertices with weights 0 or 1, all cutoff values −1, 0,
  1, and all nonnegative multiplier values 0, 1, 2), and a local harmonic
  sublevel-flux example;
- 5,832 restricted near-energy bounds and interval-integral comparisons
  (all four-gap lists with entries 1–3, exponents 3–6, thresholds 1–3,
  and cutoff radii 1–6).

The truncation test does not assume that an arbitrary finite configuration
is an equilibrium. Its sample nonzero residual is retained explicitly.
These finite checks support the proposed general argument; they do not
certify infinite rigidity, real nonintegral exponents, or the absence of
counterexamples outside the tested families. The intended positive
settlement uses an infinite analytic argument, rather than enumeration.

The implementation is independent and uses the mathematics of the cited
paper; none of the source archive's computer programs was read or copied.
