# Harmonic steps and sign consensus

Run the check with Python 3.12 or later and SymPy 1.14.0. The recorded stdout is in [check.out](check.out).

    python3 check.py

The check evaluates equations (18) and (19) on the two-node graph with both directed edges of weight one. It checks that both forms preserve the sum and that sign consensus is stationary, as in Lemma 7 of the source. The family calculation uses an exact positive-parameter substitution covering \(0<\delta<2a<1\). It checks the first iterate, both gradient signs, stationarity for an arbitrary subsequent step size and the objective gap. [arXiv:2610.24431v1](https://arxiv.org/abs/2610.24431v1)
