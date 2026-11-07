# Discrete Rogers–Shephard entropy checks

These checks concern Section VIII-B, Eq. (63) of Bang-Xian Han, Guoqing Zhan,
Chiyu Zhou and Zhuo-Nan Zhu, *Sharp Rogers--Shephard Inequalities for Rényi
Entropy*, [arXiv:2610.15716v1](https://arxiv.org/abs/2610.15716v1).

Run `python3 check.py` with Python 3 and SymPy. The expected last line is
`ALL_CHECKS_OK`; the stored output is `check.out`.

The fidelity check reproduces the source's geometric difference-law formula
and entropy increment by exact generating-function algebra. The structural
check enumerates every positive integer weight vector of lengths 1–6 with
weights in 1–4, selects the 282 log-concave vectors (160 nonmonotone), and
checks every antidiagonal central-selection identity and every majorization
prefix against the distribution with two copies of each input mass divided
by two. It independently encloses each Shannon increment below log(2).

Additional checks cover a point mass, a uniform law, a binomial law, an
asymmetric unimodal law, truncated geometric laws, and the infinite geometric
entropy formula. For parameter 999/1000 the certified gap to log(2) is between
0.000000125125 and 0.000000125126.

All evidence uses integers, exact rational numbers or symbolic algebra.
Logarithms are enclosed by a positive atanh series with an explicit rational
remainder bound; printed decimals round the endpoints outward using integer
arithmetic. No floating point calculation contributes evidence.

The finite search and formula checks support the candidate argument. They do
not establish the universal inequality, the countable-support extension, or
an equality classification. Those require a general proof.
