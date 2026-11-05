# Concave MaxTSP expectation gap

Exact supporting checks for Sergio Caracciolo and Carlo Vanoni, “Exact cost of the Euclidean Maximum Traveling Salesman Problem in one dimension”, arXiv:2610.26536v1, submitted 19 October 2026, Equation (6.30). Source: https://arxiv.org/html/2610.26536v1.

Run with Python 3, without third-party packages:

```sh
python3 check.py
```

The stored output is `output.txt`. The program uses integers and rational numbers exclusively. For costs with exponent 1/2, it compares sums of square roots using certified dyadic rational intervals derived from integer square roots. It uses no floating-point evidence.

- Fidelity: the beta-density integral reproduces the source's rank-gap moment formula (6.4), for N from 3 through 14, all ranks 1 ≤ k < N, and p in {1/5, 1/3, 1/2, 2/3, 4/5}: 450 exact cases.
- Fidelity: the edges in (6.23) form a Hamiltonian cycle with rank-gap multiplicities (m−1, 2, m−1) at (m−1, m, m+1), where N=2m. Their expected sum agrees with (6.24).
- Proposed argument: for every even N from 4 through 40 and the five rational exponents, the exact matching envelope satisfies

  N μ_m − E[E(h̄)] = μ_m p(1−p)(m−1)/(m(m−1+p)) ≤ 2p(1−p)/N ≤ 1/(2N).

  The rational-function coefficient identity is also checked by symbolic integer polynomial expansion. Here μ_k is the moment (6.4).
- Pointwise support: at p=1/2, enumerate every undirected Hamiltonian tour and every perfect matching for all ordered point subsets of the integer grids 0..6 (N=4), 0..7 (N=6), and 0..8 (N=8). The crossing matching is optimal; every tour has cost at most twice that matching cost; the reference cost is at most the best tour cost. A common positive rescaling places each configuration in [0,1] without changing these homogeneous comparisons. Totals: 24,465 tours and 1,470 matchings over 72 configurations.

The program does not integrate the maximum tour cost. It checks source fidelity, the proposed general upper-envelope algebra, and finite pointwise instances. It is not a proof of (6.30) or a Lean probe. The proposed positive settlement needs a general matching argument and order-statistic integration, not positive exhaustive enumeration.
