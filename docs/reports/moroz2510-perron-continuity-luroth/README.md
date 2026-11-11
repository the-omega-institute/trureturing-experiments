# Perron continuity: Lüroth family

Run `python3 docs/reports/moroz2510-perron-continuity-luroth/check.py`.
The program uses only integers and `fractions.Fraction`. Its output is in
`check.out`.

The fidelity check reproduces the equal-diameter cylinder formula in
Proposition 2.5 of Mykola Moroz, arXiv:2510.03373v2, for the Lüroth parameter
choice, and matches affine branch compositions to the printed series.
The remaining checks evaluate the proposed counterexample at 4/7 and the
family (7j−4)/(7j(j−1)) for 2 ≤ j ≤ 20. Periodic alternating codes stay strictly
inside their branch intervals. Exact approaching sequences and error bounds
support the proposed jump limits.

Finite rational checks do not prove the universal family or the full one-sided
limits. No Lean proof or probe is included. The stronger general classification
of the discontinuity locus is outside this report.

Source: https://arxiv.org/abs/2510.03373v2, Section 2 and the unnumbered
conjecture in Section 4.1. Related prior partial results: Min Woong Ahn,
arXiv:2610.00946, Theorem 3.5 and Remark 3.7.
