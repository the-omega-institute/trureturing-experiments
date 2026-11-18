# Markov bases for small conditional tables

Run `python3 check.py` from this directory with Python 3.12 or newer. The script uses the standard library and exact arithmetic. The recorded run is in `check.out`.

The fidelity check rebuilds the denominator row from the counts in Table 1 of Slavković, Zhu and Petrović, arXiv:1401.1397v1, §4.1. It reproduces the nine positive solutions to the source's equation. The two family checks use consecutive weights for sizes 3 through 12 and Bresinsky's four-generator weights for parameters 2 through 12. For each row they verify positive cell counts, conditional dependence of the response and omitted variables, and the exact minimum size of a Markov basis.

At a positive degree, an index is active if the residual after subtracting its weight is attainable. Active indices are adjacent if the residual after subtracting both weights is attainable. The components count the groups of factorizations linked by steps that share a positive coordinate. Each such step is a translate of a lower-degree step. An edge crossing two components has disjoint endpoint supports, so its signed difference can cross components at only one degree and for only one unordered pair of endpoints. Induction on the degree gives a minimum of one fewer move than the number of components at each nonempty fiber, summed over the degrees.

The script divides out the common gcd. A run of attainable degrees of length equal to the smallest normalized weight certifies that all later degrees are attainable. Adding the maximum pair sum to the start of this run gives a cutoff above which every pair of active indices shares a factorization. The degree scan therefore gives the exact minimum for each tested row.

Sources: [Slavković–Zhu–Petrović, Conjecture 3.3](https://arxiv.org/abs/1401.1397v1) and [Herzog–Stamate, §3.2](https://arxiv.org/abs/1308.4644v3).
