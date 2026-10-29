# Mittag-Leffler reverse modulus inequality

GGPS Conjecture 25 in arXiv:2410.11852v4 asserts, for 1 < α < 2 and β ≥ 0,
that the global inequality |E_{α,β}(z)| ≥ E_{α,β}(Re z) is equivalent to
α−1 ≤ β ≤ α, and also to complete monotonicity of the reciprocal on (0,∞).
Here E_{α,β}(z) = Σ_{n≥0} zⁿ/Γ(αn+β).

Run `python3 check.py`. Python 3's standard library suffices. The stored output
is `check.out`; a successful run ends with `ALL_CHECKS_OK`.

The program uses integers and `fractions.Fraction` throughout. It certifies
strict counterexamples with α = 3/2, β = 1/2, and z = r(−1+i/2), for each
r ∈ {1024,1536,2048}. In each case it encloses the actual defining infinite
series, proves E_{3/2,1/2}(−r) > 0, and compares a rational upper bound for
|E_{3/2,1/2}(z)|² with a rational lower bound for E_{3/2,1/2}(−r)².
These are finite exact certificates of failure of the source assertion.

The fidelity check reproduces the source's formulas E_{1,1}(z)=exp(z) and
E_{2,2}(z)=sinh(√z)/√z at the level of coefficients n=0,…,80. The second
formula is read at z=0 by its removable value. This finite coefficient check
is distinct from a proof of the analytic identities.

## Certificate method

For s a positive integer or half-integer, 1/Γ(s) has the exact form a+b/√π
with a,b nonnegative rationals. The half-integer formula is
1/Γ(k+1/2) = 4ᵏ k!/((2k)!√π). Machin's identity
π=16 arctan(1/5)−4 arctan(1/239), alternating-series bounds with 160 terms,
and integer square roots give rational bounds for 1/√π.

After terms n=0,…,1199, the two parity tails are bounded independently.
For s=(3n+1)/2, the quotient between terms n and n+2 is
z²/[s(s+1)(s+2)]. Put R=|Re z|+|Im z|. The absolute quotient is at most
q=R²/[s(s+1)(s+2)], and decreases along each parity. Once q<1, the tail
starting at n is bounded by (a+b)Rⁿ/(1−q), using 1/√π<1. The sum of both
parity bounds is less than 10⁻⁶⁰ at each displayed real and complex input.
The reported rectangles are rounded outwards to a rational grid. There is
no floating-point special-function evaluation or unbounded-tail omission.

## Analytic family and limits of the check

The sector expansion in Haubold–Mathai–Saxena, arXiv:0909.0230v2, §6,
(6.9) and (6.11), gives at β=α−1

E_{α,α−1}(z) = c_α z⁻² + O(|z|⁻³),

where c_α=−1/Γ(−α−1)=α(α−1)(α+1)/Γ(2−α)>0. The linear coefficient
vanishes because 1/Γ(−1)=0. Along a ray z=r(−1+it) outside the exponential
sector, the quotient of the complex modulus and the positive real-axis
value tends to 1/(1+t²)<1. For each fixed α in (1,2), choose such a ray,
a sufficiently large finite r, and then use continuity in β at that fixed
complex input to obtain an interval [α−1,α−1+ε_α] of violations, with
0<ε_α<1. The same input works throughout that interval.

The program independently checks the normalized coefficient identity
α(α−1)(α+1)=u(u+1)(u+2), u=α−1, and the exact ray factor 4/5 for t=1/2.
It does not verify the sector expansion, the all-α limit argument, parameter
continuity, or an explicit interior ε_α. Those analytic steps belong to the
preregistration's ordinary mathematical outline. No Lean proof or kernel
certification of the new conclusion is supplied by this experiment.

Sources: Garrappa–Gerhold–Popolizio–Simon, *On some inequalities for the
two-parameter Mittag-Leffler function in the complex plane*, JMAA 551(1)
(2025), 129588, DOI 10.1016/j.jmaa.2025.129588, arXiv:2410.11852v4;
Haubold–Mathai–Saxena, *Mittag-Leffler Functions and Their Applications*,
arXiv:0909.0230v2. The check program is an independent implementation.
