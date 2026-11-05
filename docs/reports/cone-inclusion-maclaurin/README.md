# Cone inclusion and Maclaurin means

Exact bounded checks for the unnumbered Conjecture after equation (1.5) in
Ning Cao, Xuanxuan Hu and Feida Jiang, *New Newton-type inequality and a unified
perspective on Maclaurin's inequality via the inclusion relations of cones*,
[arXiv:2610.28544v1](https://arxiv.org/abs/2610.28544v1).

Run with Python 3, using only its standard library:

```sh
python3 docs/reports/cone-inclusion-maclaurin/check.py
```

`check.out` is the output. Every arithmetic operation carrying evidence uses
integers or `fractions.Fraction`. Radical comparisons raise nonnegative rational
bases to a common positive integer power; the program never evaluates roots or
uses floating point.

The fidelity check reproduces the inequality in Theorem 1.1 at n=3, k=2,
x=(1,2,3), and also tests its inequality on the specified finite grid. The endpoint
identities Pi_1=sigma_n and Pi_n=sigma_1 are checked independently.

The settlement checks verify all proposed obstruction vectors for n=2 through
12, all indices, and the one-negative-coordinate coefficient formula at five
rational parameter values in each dimension. Elementary symmetric functions are
computed by an independent polynomial-coefficient recurrence. The four proposed
index conditions are tested against all 1,708 coordinate multisets of dimensions
2 through 6 from {-2,-1,-1/2,0,1/2,1,2}. All 18,764 applicable ordered cone/mean
comparisons pass, including zero and negative coordinates whenever admissible.

The program tests Gamma's algebraic closed domain sigma_j>=0, and P's closed
domain of nonnegative coordinate sums. Identification of the former with the
topological closure belongs to the mathematical proof. These finite checks do
not prove the universal conjecture, continuity, or the inclusion classification;
those remain obligations for a later Lean delivery. The expected proof uses the
source's Corollaries 1.2 and 1.4, Theorem 1.3, and classical Maclaurin; it must
formalize those ingredients rather than treating them as axioms.
