# Sparse two-round one-excellent-element search

Exact finite checks for Nagananda K G and Jong Sung Kim, *Statistical two-round search for one excellent element*, arXiv:2605.15612v1, Section 4(i), with the model of Section 1.1 and Eqs. (1)–(2). The proposed general result is

\[
\lim_{n\to\infty} M_n^{(2)}(\lambda,\alpha)/\ln n
= (1-\alpha)/\ln 2,
\qquad \lambda>0,\quad e^{-\lambda}<\alpha<1.
\]

Run with Python 3 and its standard library:

```sh
python3 docs/reports/nagananda-kim-2026-two-round-search-limit/check.py
```

`output.txt` contains the complete output. All arithmetic is integer or `fractions.Fraction` arithmetic. Logarithmic integer thresholds use bit lengths and integer square roots.

The fidelity checks reproduce the source's binomial excellence law, finite-population empty-event feasibility bound, and expected-cost identity for a global existence test followed by a second batch (using singleton tests here). They cover 24 parameter pairs.

The dynamic programme enumerates every reduced batch for n = 1, 2, 3, and every possible first-response branch. Removing empty or repeated tests and changing their order cannot improve information. The decoder chooses the index with largest unconditional mass of excellence in each final response cell. Reporting an index without testing it individually is allowed by the source. Pareto pruning preserves the exact deterministic optimum for every success threshold. Taking pairwise rational mixtures gives the randomized convex envelope; arbitrary independent randomisation cannot improve on that finite convex hull. The output checks 33 thresholds against the proposed exact tail-sum converse. For example, with n = 2, p = 1/4 and required success 21/64, the deterministic optimum is 1 and the randomized optimum is 5/12. The distinction is material at finite n.

The proposed block-and-binary-label construction is checked by enumerating all excellence configurations for every n from 2 through 14, at two rational parameter pairs. All 26 cases meet the unconditional success constraint and the exact cost and collision bounds. The balanced-block pair-count bound is also checked on 4,160 integer cases.

These finite checks support the planned proof ingredients. They do not establish an asymptotic limit. The settling argument is a general tail-sum lower bound and a construction with a vanishing collision loss, rather than a positive finite enumeration.

The source's probability parameter is defined only for n > 0 and 0 ≤ λ/n ≤ 1. Every checked parameter satisfies this condition. The bound construction is asserted only when its success threshold is attainable; infeasibility would be reported by returning no construction.

Source: https://arxiv.org/abs/2605.15612v1. The code is independently written; no paper text or third-party code is reproduced in this report directory.
