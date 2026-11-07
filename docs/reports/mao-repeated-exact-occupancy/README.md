# Repeated exact-occupancy MAO inequality

`python3 check.py > output.txt` uses Python 3's standard library, integer counts
and `fractions.Fraction`. No floating point is used as evidence.

The source is Mao–Xue, arXiv:2512.06880v1, Theorem 2 and appendix Conclusion,
and arXiv:2608.23998v1, Definitions 2–4, Theorems 5–6, and the section 7
Remark after equations (92)–(93).

Fidelity: the script reproduces the section 8 exhaustive example with population
size 9 and sampled subset sizes (3,2,2,1). It enumerates all 979,776 original
configurations and independently computes the source's transversal sums.
The first four factorial and raw moments and the variances agree exactly for
exact and threshold occupancy levels 1–4. The fully occupied inequality is also
checked for all orders 1–9 for this model.

Refutation candidate: population size and moment order 2, four independent
uniform singleton subsets, and exact occupancy 2 give both the mean and the
second factorial moment equal to 3/4. Thus 3/4 > (3/4)^2 = 9/16. The
variance is 15/16. All subset sizes are positive and proper.

For 2k independent uniform singleton subsets of a two-element population,
the exact-half indicators coincide. The outlined family has common mean and
second factorial moment u_k = binom(2k,k)/2^(2k−1). The program checks all
central ordered pattern pairs through k=5, all original samples through k=7,
and the exact recurrence u_(k+1)(2k+2) = u_k(2k+1) through k=100.
The endpoint k=1 gives equality; all tested k≥2 give strict reversal.

These are finite exact checks. The unbounded complement-bijection and
induction arguments require separate proof; no Lean proof is supplied here.
