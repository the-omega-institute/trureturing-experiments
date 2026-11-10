# Strong stochastic dominance in the unequal-probability coupon collector

The target is the unnumbered conjecture following Definition 8.2 in Zoroa,
Lesigne, Fernández-Sáez, Zoroa and Casas, arXiv:1608.00155v1, retained after
Definition 7.2 of J. R. Soc. Interface 14 (2017), 20160643,
DOI: 10.1098/rsif.2016.0643. This is a preregistration experiment, not a Lean
proof or a claim of a kernel-verified settlement.

Run `python3 docs/reports/1608.00155-strong-dominance/check.py` from the
experiments repository. The program needs only Python's standard library.
Its stored stdout is `check.out`; successful execution ends in `ALL_CHECKS_OK`.
All computations use integers or `fractions.Fraction`.

The fidelity checks reproduce the source's inclusion-exclusion PMF and mean
identities in 600 cases (1 ≤ N ≤ 4, denominator 6, 1 ≤ n ≤ 5, including zero
weights). They also verify the source's stated strong-order result on 12,012
majorized ordered pairs for N = 4, n = 2, 3, 4 and denominator 6. These finite
checks do not establish the source's full universal special-case results.

Independent enumeration of all 4^5 = 1024 words gives a counterexample at
N = 4, n = 5: p = (48,48,2,2)/100, q = (50,46,3,1)/100. With denominator
100^5, the PMF difference at k = 1,2,3,4 is
(-8855220,1524660,-6565440,13896000). Its positive entry at k = 2 followed
by a negative entry at k = 3 excludes the source's threshold convention
(k ≤ c on the low side and k > c on the high side). Proper cumulative
probability differences are negative, consistent with ordinary stochastic order.

Exact polynomial enumeration independently checks the singleton/pair/product
counting formulas for p_t = (1/2-2t,1/2-2t,2t,2t) and
q_t = (1/2,1/2-4t,3t,t). Coefficientwise rational interval bounds certify
positive normalized sorted vectors, p_t ≺ q_t and PMF sign pattern (-,+,-,+)
for every real 0 < t ≤ 1/100. The output includes all four difference
polynomials and the exact enclosure of each polynomial after removing t^2.
The t = 0 boundary is excluded; there p_t = q_t.

The source allows nonnegative probabilities. The witness and the whole stated
interval have strictly positive probabilities, so the obstruction persists
under the positive-probability restriction as well. The proposed Lean delivery
uses the concrete instance to refute the universal claim; the interval explains
the mechanism and is not a requirement to add an extra Lean theorem.
