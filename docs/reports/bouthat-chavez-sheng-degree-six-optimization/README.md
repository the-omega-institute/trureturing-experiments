# Degree-six iid graph-norm optimization

The proposed answer to Bouthat, Chávez and Sheng, *A general framework for inequalities on simple graphs*, arXiv:2606.17356v1, Section 5(ii), reduces both extremal norm ratios to two explicit one-parameter families. For standardized laws write t = E(Y^4), a = E(Y^3)^2 and u = E(Y^6). The squared-polynomial certificates give 0 ≤ a ≤ t − 1 and z = u − t² − a ≥ 0. For every trace-zero spectrum its sixth moment decomposes as (1 − θ)S(t) + θT(t) + z p₆, where θ = a/(t − 1) for t > 1. Both boundary families are actual finite-support laws. The same coefficients apply to numerator and denominator, giving a convex combination of ratios. Two common endpoints and at most four stationary points determine the extrema.

Run with Python 3 and SymPy:

```sh
python3 docs/reports/bouthat-chavez-sheng-degree-six-optimization/check.py
```

The stored output is `check.out`. All mathematical calculations use integer, rational or exact symbolic arithmetic. The checks comprise:

- 200 exact evaluations of the source's Rademacher comparator formulas, including the sixth-moment formula;
- independent convolution of iid atomic sums for all 24 nonzero, trace-zero integer multisets of length 2–5 with entries −2 through 2, against 51 probability laws;
- the decomposition and ratio bounds for every one of the 1094 labeled nonempty-edge graphs of orders 2–5, against matching, complete, one-edge and K₄ comparators: 223176 graph/comparator/law combinations;
- symbolic squared-polynomial identities, boundary-family moments, the conical decomposition and derivative formula;
- the interior optimizer for K₁,₃ versus K₄ at t = (7 + 2√511)/19, strictly exceeding the Rademacher and Schatten-limit values.

These are exact finite checks and symbolic identities, not a Lean proof of the universal optimization theorem. The source's definition is the literal expectation divided by Γ(7) = 720. The malformed left-hand side of its first Section 3.4 display is not used as a norm definition, and its displayed bipartite maximum is not used as an unpowered norm value. The fidelity checks use the unambiguous Rademacher formulas in Section 3.5.
