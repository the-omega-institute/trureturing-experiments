# Nonelliptical GFT flatness: exact check

`check.py` checks the covariance functional in Sections 3.2 and 5.1 of Argyn
Kuketayev, *The Sampling Distribution of the Log-Euclidean Distance Between
Sample Correlation Matrices*, arXiv:2608.25288v1, Sections 5.2–5.3 and the final
question in Section 6.5. The source is available at
https://arxiv.org/abs/2608.25288v1. No source text or source figures are redistributed.

Run with Python 3 and SymPy 1.14.0:

```sh
python3 docs/reports/kuketayev-log-euclidean-flatness/check.py
```

`output.txt` is the exact output. The script uses integer, Fraction and symbolic
exact arithmetic throughout. It reproduces the source's bivariate Gaussian
identity V(r)=1 from Wick fourth moments. It then independently expands the
source influence function, verifies a 72-label atomic law's moments for every
t>0, and obtains v_t(0)=1 and v_t'(0)=6/(t²+5t+4)>0. At t=1/2 the derivative is
8/9. A rational parametrization rho(h)=h/(1+h²/4), |h|<2, yields a rational
variance with the same derivative. Independent symmetric ±1 coordinates give
an identity baseline covariance in every dimension; complete covariance matrices
are checked for dimensions 2 through 8. Symbolic resolvent block isolation is
checked in dimensions 2, 3 and 4. The small dimension checks do not replace the
general block argument.

The finite noise takes 0 and ±sqrt(3t), with probabilities 2/3 and 1/6 each,
so its first four moments match the Gaussian noise of the candidate construction.
The law has finite support and satisfies the source's finite-fourth-moment
assumptions; it makes no smooth-density claim. Its atom locations depend smoothly
on the local parameter. At t=1/2 all latent moments are rational; after rational
parametrization all moments at rational h are rational as well.

The result is an exact algebra check supporting a preregistration outline. No
Lean proof, statistical simulation, or claim of exhaustive literature novelty
is part of this artifact.

Code and report: Apache-2.0, under the repository license. Source-paper rights remain
with its author and the arXiv distribution license.
