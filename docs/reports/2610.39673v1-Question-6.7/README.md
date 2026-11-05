# Random semigroups at critical density

Source: Guy Blachar, Be’eri Greenfeld and Margaret Shideler, *Random semigroups á la Gromov*, arXiv:2610.39673v1, Question 6.7 (29 October 2026).

Run `python3 check.py`. The program uses integers and rational arithmetic throughout. It checks 16 instances of the source's ball/pair counts, 120 exact endpoint-pattern probabilities against complete enumeration (including Proposition 6.6's avoidance formula), and 63 evaluations of the proposed factorial-language mass formula against all 16,384 marking configurations on binary words of length at most three. The latter checks only closed languages with maximum word length at most two.

For n=3 it encloses the proposed limiting triviality probability `(1-exp(-4/3))^3` between `399343300/1000000000` and `399343301/1000000000`. The exponential enclosure uses a positive Taylor series and a rational geometric tail bound.

The seeded simulation samples exactly 3^l distinct uniformly chosen unordered relations for each of 160 trials at each l from 2 through 6. Congruence closure on words of length at most l+1 adds only equalities derived from the original relations. Finite classifications require a complete multiplication table, associativity, and evaluation of all relations. A failed finite certificate is unresolved, and is not classified as an infinite or nontrivial semigroup. Finite certificates are algorithmic checks, not Lean proofs.

| l | finite certified / trials | trivial certified / trials | unresolved / trials |
|---|---|---|---|
| 2 | 151/160 | 89/160 | 9/160 |
| 3 | 156/160 | 56/160 | 4/160 |
| 4 | 156/160 | 68/160 | 4/160 |
| 5 | 158/160 | 67/160 | 2/160 |
| 6 | 160/160 | 56/160 | 0/160 |

The output includes the full empirical size histograms and exact finite-l probabilities that every generator occurs as a relation endpoint. Triviality and endpoint marking disagree on two certified samples at l=2 and one at l=3; no disagreement occurs in the certified samples at l=4,5,6. The proposed equivalence is asymptotic, subject to the bounded-nilpotency and no-short-to-short-relation reduction; it is not asserted for every finite presentation.

The deterministic monomial reduction is checked on all 694 simulated samples satisfying its hypotheses, including every product of surviving words. These finite experiments do not establish convergence or any unbounded theorem. The proposed proof must supply finite-dimensional asymptotic independence, the uniform reduction using source Theorem 6.2, and the probability-one finiteness and normalization of the limiting language law.
