# Exact finite checks for mixed multiplicative convergence

Run `python3 check.py` with Python 3 and SymPy 1.14.0. The stored output is `check.out`.

The script reproduces the differential-operator identity in Section 3.1 of Katsunori Fujie, *Regularity and Convergence Properties of Finite Free Convolutions*, arXiv:2505.15575v1. It checks the permutation formula and finite atom lower bounds on every multiset pair of degrees one through four over the displayed alphabets. It checks the one-root replacement bound on seeded rational examples of every degree one through twelve, in both input positions, including zeros and repeated roots. Root counts use exact rational Sturm isolation; no floating-point roots or tolerance decisions are used.

The final zero-multiplicity example explains why only a lower bound is used when the first polynomial has signed roots. These finite checks provide no proof of either asymptotic convergence assertion and no claim about priority or journal-version question status.

Source: https://arxiv.org/abs/2505.15575v1 .
