# Rainio–Kargar logarithmic-mean metric check

Source: Oona Rainio and Rahim Kargar, *Intrinsic metrics defined with arithmetic
and logarithmic mean values*, [arXiv:2308.06576v1](https://arxiv.org/abs/2308.06576v1),
Conjecture 4.6 and equations (4.1)–(4.2).

Run `python3 check.py` with Python 3 and SymPy 1.14.0 installed. The program uses
integers, rational numbers, exact radicals, and symbolic identities. It performs
no floating point evaluation. Successful output is stored in `check.out` and ends
with `ALL_CHECKS_OK`.

The fidelity check reproduces the unbounded-function assertion of source Lemma
4.5: at the ball points ±e₁/2 and 0, the unbounded triangle defect is
2(1 − log 2)/c ≥ 1/(2c) > 0 for every c > 0.

For the bounded transform in Conjecture 4.6, the same points at c = 2 have defect
at least 1/57. The symbolic family uses r = (c − 1)/(2c − 1) for every c > 1.
The half-space family uses (±3e₁ + eₙ)/4 and eₙ for every c ≥ 4. The program checks
the source substitutions, exact Euclidean lengths, positive rational lower bounds,
and strict radical comparisons. Zero-padding checks use dimensions 2, 3, 4 and 7;
the written dimension-independent construction applies to every n ≥ 2.

The logarithm estimate is certified independently by the identity between the
chord of 1/t and 1/t: their difference is (t − a)(b − t)/(abt) on [a,b]. Its
nonnegativity and exact integral give the trapezoidal bound used here. This is
the classical inequality already supplied by the repository's
`D5.S3.Divergence.MeanKernels.LogarithmicMeanSandwich.logMean_kernel_sandwich`.
The program checks the algebraic and derivative/integral certificates; it is not
a Lean proof or an independent formal verification of real analysis or topology.
The identification of the actual topological boundary distances, the source
binding, and the final negation of the metric claim remain for Lean verification.

These results support an expected refutation. They make no assertion about
metricity on the ball for 0 < c ≤ 1 or on the half-space for 0 < c < 4.
