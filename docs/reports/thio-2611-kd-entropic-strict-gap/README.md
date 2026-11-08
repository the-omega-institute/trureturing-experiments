# KD entropic strict gap: independent finite checks

The mathematical target is Appendix B, equations (63) and (67), of Thio, De Bièvre, Barnes and Arvidsson-Shukur, arXiv:2611.04009v1. These files check the proposed tangent-descent route; they do not prove the probability-one assertion.

Run with Python 3.12 or newer and SymPy 1.14:

```sh
python3 check.py
```

`check.py` uses only exact rational and symbolic arithmetic for its certificates. Natural logarithms of positive rational numbers are enclosed by the convergent atanh series with an explicit geometric-tail bound. Zero entropy terms are zero. The checks reproduce the source's three-dimensional KD matrix and marginals, all basis-state KD laws for its matrix (58), and the displayed four-dimensional KD-positive state for U_0.

Nine deterministic pseudorandom complex rational Cayley unitaries, three each in dimensions 2, 3 and 4, satisfy the exact KD-real nullity 2d−1. Every overlap is nonzero. For each, the program certifies the least row/column entropy, nonzero logarithmic variance, the symbolic derivative identity −2 Var(log p), and a normalized rational pure-state witness with entropy sum below that least vertex entropy by more than 2^−60. Exact entropy enclosures also put the B finite-difference error relative to −2 Var(log p) below 1/1000 and A(t)/t below 1/1000. The tangent used for the finite witnesses approximates each logarithm within 2^−63; it is stated as an approximation, while the derivative identity uses symbolic logarithms. The KD-real nullity is the deterministic geometric criterion from Langrenez et al., arXiv:2405.17557v2, connecting these examples to the basis-projector polytope.

The stored `check.out` additionally includes the optional diagnostic:

```sh
python3 check.py --haar-diagnostics
```

That option needs NumPy. `haar_diagnostics.py` generates nine floating-point complex Gaussian QR Haar samples in dimensions 2, 3 and 4, searches for feasible pure-state entropy upper bounds by local optimization, and compares finite differences with −2 Var(log p). Those lines are explicitly labeled FLOAT and are exploratory diagnostics only. They are not evidence for any exact inequality, any global optimum, the KD-positive polytope classification for those samples, or a Haar-almost-sure statement. The exact certificates above use no floating-point arithmetic. Both modes end with ALL_CHECKS_OK on success.

No Lean proof or kernel settlement is included. The all-dimensional argument, the generic KD-polytope input and the Haar-null exceptional-set argument remain formalization obligations.
