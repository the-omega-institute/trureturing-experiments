# A prime-weighted continuation using RH

Let a finite family have distinct odd moduli greater than one and arbitrary
fixed residues. Write the distinct prime factors of its LCM in increasing
order. Each row gives a sufficient condition for an uncovered integer.
All prime-power heights, including those of the head primes inside mixed
moduli, are unrestricted. The number of prime factors above the indicated
cutoff is unrestricted.

| Cutoff \(B\) | Maximum number of prime factors at most \(B\) | Additional condition |
| ---: | ---: | --- |
| 449 | 7 | Every such family |
| 2543 | 8 | Every such family |
| 8231 | 9 | If the ninth prime is at most \(B\), it is at least 271 |

These are ordinary proofs using the published RH result and the existing
survivor measure constructions. Alpöge's Lean declaration has type
\(\mathrm{RiemannHypothesis}\) without hypotheses. The independent
statement comparison and axiom audit verify that declaration against
Mathlib. The analytic step below only needs ordinary RH.
[Proof declaration](https://github.com/the-omega-institute/riemann-hypothesis/blob/1ef6c4e2e944fb2784ecf3168a84810bf50ca6db/src/Main.lean) ·
[Statement comparison](https://math.mit.edu/~drew/grh/con-leche-compare.log) ·
[Axiom audit](https://math.mit.edu/~drew/grh/audit.log)

### Survivor measures and continuation

Take the actual primes at most \(B\) as the head. Resolve their finite
coordinates at the full heights required by the original family. The
head submeasure avoids the classes supported entirely in the head. Its
mass is at least \(m\), and the integral of the fourth power of every
complete divisor query is at most \(K\). Such a query counts hits against
one fixed residue for each divisor of the head period, including the
unit. Hölder bounds the integral of any product of four queries by the
same \(K\) on this measure. The following source bounds apply, with
unused head coordinates padded when necessary.

| Head | \(m\) | \(K\) | Source |
| --- | ---: | ---: | --- |
| At most seven primes | \(7235955529/450000000000\) | \(82461631627375/127401984\) | [Seven-prime source and fourth moment, HM1–HM6](https://github.com/the-omega-institute/trureturing-experiments/blob/80d407671ff57ec39cec5e7244b2fee1df2059c3/docs/reports/erdos7-odd-covering/profile-notes/arithmetic/700-749/734-seven-and-eight-full-height-heads-admit-quartic-prime-tails.md) |
| At most eight primes | \(1/33750\) | \(590500\) | [Capped eight-prime source, (1)–(5)](https://github.com/the-omega-institute/trureturing-experiments/blob/80d407671ff57ec39cec5e7244b2fee1df2059c3/docs/reports/erdos7-odd-covering/profile-notes/arithmetic/750-799/779-the-original-capped-source-lowers-the-general-eight-prime-tail-cutoff.md) |
| Nine primes, ninth at least 271 | \(1/2500000\) | \(450000\) | [Nine-prime source, (1)–(5)](https://github.com/the-omega-institute/trureturing-experiments/blob/80d407671ff57ec39cec5e7244b2fee1df2059c3/docs/reports/erdos7-odd-covering/profile-notes/arithmetic/750-799/780-a-ninth-prime-from271-admits-an-unrestricted-larger-tail.md) |

For the last row with at most eight head primes, use the eight-prime
source and the cutoff of 2543. In the nine-prime case the separation
condition is imposed on the actual ninth prime.

Process the larger actual primes in increasing order and charge each
remaining class at its largest prime. For \(t=1/(p-1)\), put
\[
 A(t)=15t+50t^2+60t^3+24t^4,\quad
 c_\delta=\frac{27}{256\delta^3(1-\delta)},\quad
 g_\delta(t)=1+\frac{A(t)}{1-\delta}.
\]
The fourth-moment continuation bounds the absolute loss at a prime by
\(Kc_\delta t^4\), with next fourth-moment allowance \(Kg_\delta(t)\).
It uses the actual forbidden fraction \(\beta(x)\) on the \(p\)-fibre
over each old point \(x\), restricting Haar to the allowed points and
dividing by \(1-\min(\beta(x),\delta)\). The mass loss is thus the integral
of \((\beta-\delta)_+/(1-\delta)\). All old cofactors and phases remain
those of the original classes. At each new exponent, distinctness gives
at most one phase per old divisor. This is the hypothesis used in the
complete-query estimate. The inequalities hold for each \(0<\delta<1\),
so the clipping parameter can be fixed separately for each prime before
processing the family. [Quartic continuation, HM7–HM13](https://github.com/the-omega-institute/trureturing-experiments/blob/80d407671ff57ec39cec5e7244b2fee1df2059c3/docs/reports/erdos7-odd-covering/profile-notes/arithmetic/700-749/734-seven-and-eight-full-height-heads-admit-quartic-prime-tails.md)

### Every prime above the computation range

Set \(T=200000\) and use \(\delta=2/5\) for \(p>T\). The binomial
coefficients give \(g_{2/5}(t)\le(1+t)^{25}\) for \(t\ge0\).
For \(x\ge T\), the reciprocal prime-product bound gives
\[
 \prod_{T<p\le x}\frac p{p-1}
 \le \frac{243}{241}\frac{\log x}{\log T},
 \qquad 11<\log T<13.
\]
Indeed \(3^{11}<T\) and \((27/10)^{13}>T\). The constant follows by dividing
the upper bound at \(x\) by the lower bound at \(T\) and using
\(\log T>11\). [Rosser–Schoenfeld, Theorem 8](https://github.com/the-omega-institute/trureturing/blob/44bb5611a7f789fedcf93c1d53a251ab39837f1a/Library/Arith/rosser1962approximate.md)

With \(\vartheta(x)=\sum_{p\le x}\log p\), write
\(E(x)=\vartheta(x)-x\). Under RH,
\(|E(x)|\le\sqrt{x}(\log x)^2/(8\pi)\) for \(x\ge T\), by weakening
Lee–Nosal's explicit bound. Since \((\log x)^2/\sqrt{x}\) decreases
for \(x\ge T\), \(|E(x)|\le\varepsilon x\) throughout this range, where
\(\varepsilon=169/(24\cdot447)\). [Lee–Nosal, Theorem 1.2(3)](https://arxiv.org/abs/2312.05628)

Let \(f(x)=(\log x)^{24}/x^4\). It is positive and decreasing over
\([T,\infty)\). Partial summation over the primes strictly greater than \(T\)
gives
\[
 \sum_{p>T}\log(p) f(p)
 = \int_T^\infty f(x)\,dx-E(T)f(T)-\int_T^\infty E(x)f'(x)\,dx
 \le (1+\varepsilon)\int_T^\infty f(x)\,dx+2\varepsilon T f(T).
\]
The endpoint at infinity vanishes since \(xf(x)\to0\).
In the prime loss bound use
\((p-1)^{-4}\le (T/(T-1))^4p^{-4}\), then apply the product
bound to all moment multipliers from \(T\) to \(p\).
After dividing out the moment allowance at \(T\), the resulting prime
sum is at most
\[
 R=\frac{5625}{2048}\left(\frac{243}{241}\right)^{25}
       \frac{T}{(T-1)^4}
 \left[
 \frac{1+\varepsilon}{33}
   \sum_{j=0}^{25}\frac{25!}{(25-j)!\,33^j}
   +\frac{2\varepsilon}{11}
 \right].
\]
For the integral, use \(f(x)\le (\log x)^{25}/(11x^4)\) and integrate by
parts repeatedly. The integral of \((\log x)^{25}/x^4\) is
\((\log T)^{25}/(3T^3)\) times the displayed factorial sum with
\(3\log T\) in place of 33. The endpoint term uses
\(T f(T)=(\log T)^{25}/(T^3\log T)\).
These estimates bound the full infinite comparison sum and therefore
every finite choice of original primes above \(T\).

### Exact cutoffs

Put \(S=10^{24}\) and start with \(u=\lceil SR\rceil=44446336\).
Visit every prime in \((B,T]\) in decreasing order. At each such prime
replace \(u\) by
\[
 \min_{\delta\in\{1/4,\,3/10,\,1/3,\,3/8,\,2/5\}}
 \left\lceil Sc_\delta/(p-1)^4+g_\delta(1/(p-1))u\right\rceil .
\]
A minimizing parameter at each step fixes the forward clipping schedule.
The recurrence is an upper bound on the remaining loss per unit of
moment allowance, multiplied by \(S\). Missing any scheduled prime only
decreases this allowance because the loss term is nonnegative and
\(g_\delta\ge1\). All inputs and rounding operations are rational. The
table gives the integer replay values and strict lower bounds for the
final live mass \(m-Ku/S\).

| \(B\) | Primes in \((B,T]\) | Final \(u\) | Final live mass greater than |
| ---: | ---: | ---: | ---: |
| 449 | 17897 | 24477013422083484 | \(1/5000\) |
| 2543 | 17612 | 49669053146940 | \(3/10^7\) |
| 8231 | 16952 | 888884725255 | \(1/10^{12}\) |

Positive live mass gives a point on the finite CRT carrier avoiding
every original class, hence an uncovered integer. In particular, any odd
distinct cover would need at least eight support primes at most 449 and
at least nine at most 2543. If its ninth prime is at least 271, it would
also need a tenth prime at most 8231. Extending these bounds to all
possible small-prime supports requires a stronger survivor estimate
for the head.

### The second-moment threshold at 73

Suppose a probability on the full-height survivors of the head at 73 has
\(\max_L\mathbb E L^2\le139.1\), with the maximum over complete divisor
queries as above. Then the whole original family is noncovering. The
proof tracks the mass lost and the second moment on the same
submeasure. In this scalar continuation, even optimal clipping at each
prime can certify a threshold strictly below \(139.26\) when every
prime above 73 may occur. The known star family at 73 forces a second
moment greater than \(139.59\) for every supported probability.
[Star family](https://github.com/the-omega-institute/trureturing-experiments/blob/80d407671ff57ec39cec5e7244b2fee1df2059c3/docs/reports/erdos7-odd-covering/problem-details/04-a-complete-star-family-refutes-the-unrestricted-gamma-73-bound.md)

At an intermediate step write \(m\) for the live mass, \(K\) for a bound
on every complete-query second moment, and \(D\) for the actual mass
deleted by the next prime \(p\). Put
\[
 a_p=\frac{3p-1}{(p-1)^2},\qquad b_p=\frac{1}{4(p-1)^2}.
\]
The same restricted kernel gives \(D\le Kb_p/(\delta(1-\delta))\).
This follows from \((y-\delta)_+\le y^2/(4\delta)\), by completing each
old cofactor query and bounding each product by Cauchy–Schwarz. In the
enlarged second moment, the term with both new exponents zero loses at
least \(D\), since the corresponding old complete query contains the
unit term and is everywhere at least one. The other exponent pairs use
the cylinder cap \(p^{-j}/(1-\delta)\) at their maximum depth \(j\).
There are \(2j+1\) such pairs. Thus
\[
 m'=m-D,\qquad K'\le\left(1+\frac{a_p}{1-\delta}\right)K-D.
\]
The loss credit uses the same \(D\) in both bounds. [Second-moment continuation, T1–T6](https://github.com/the-omega-institute/trureturing-experiments/blob/80d407671ff57ec39cec5e7244b2fee1df2059c3/docs/reports/erdos7-odd-covering/problem-details/08-arbitrary-head-transfer-by-the-joint-load-invariant.md)

For a continuation requiring \(m'>vK'\), with \(0\le v<1\), the
preceding step can use the sufficient threshold
\[
 \Phi_p(v)=v+\min_{0<\delta\le1/2}
 \left\{\frac{b_p(1-v)}{\delta}
       +\frac{b_p(1-v)+a_pv}{1-\delta}\right\}
 =v+\left(\sqrt{b_p(1-v)}+\sqrt{b_p(1-v)+a_pv}\right)^2 .
\]
Indeed the worst loss permitted by the scalar bounds is the largest
one, since it reduces \(m'-vK'\) by \((1-v)D\). If that loss
can exhaust the live mass, the scalar bounds already allow zero
survival. The minimizer is
\(\sqrt{b_p(1-v)}/(\sqrt{b_p(1-v)}+\sqrt{b_p(1-v)+a_pv})\).
The formula increases on \(0\le v\le1/4\). The sum of the two
radicands increases there since \(a_p>2b_p\), and their product
has derivative at least \(b_p(a_p-3b_p)/2>0\).
Also \(\Phi_p(v)\ge v\) on \(0\le v<1\), so a missing prime can
be skipped. A scalar threshold that reaches \(1/4\) already
excludes any head with \(m/K\le1/4\).

For primes beyond \(T=5000000\), take \(\delta=1/2\) and use
\(1+2a_p\le(p/(p-1))^6\). Dropping the loss credit gives a
sufficient terminal bound. In the preceding partial summation use
\(f(x)=(\log x)^5/x^2\), \(15<\log T<16\) and
\(\varepsilon=16^2/(24\cdot2236)\). The log bounds follow from
\((11/4)^{15}<T<(27/10)^{16}\) and \(27/10<e<11/4\).
The prime-product and Chebyshev bounds give this terminal threshold
per unit of second moment:
\[
 R_2=\left(\frac{451}{449}\right)^6\frac{T}{(T-1)^2}
 \frac{1}{15}\left[
 (1+\varepsilon)\sum_{j=0}^{5}\frac{5!}{(5-j)!\,15^j}
 +2\varepsilon\right].
\]
The factorial sum comes from integrating \((\log x)^5/x^2\).
The boundary term again accounts for the possible negative error
in \(\vartheta(T)-T\).

Use \(S=10^{24}\) and visit all 348492 primes in \((73,T]\) in
decreasing order. For an integer \(u\), put \(X=S-u\),
\(Y=X+4(3p-1)u\), and \(d=4(p-1)^2\). Then the floor of
\(S\Phi_p(u/S)\) is
\[
 u+\left\lfloor\frac{X+Y+\lfloor\sqrt{4XY}\rfloor}{d}\right\rfloor .
\]
One more than the floor gives a strict upper rounding. The replay
starting from \(\lceil SR_2\rceil=19978401471981184\) and rounding
up at each step ends at \(7187192385814385692968\). Its ratio to
\(S\) is less than \(1/139.1\), so the stated head moment suffices.
Starting from zero and rounding down gives \(7181063643367434956708\),
whose ratio to \(S\) is greater than \(1/139.26\). This replay
bounds the best scalar threshold from below even for a tail that
ends at \(T\), with every prime in \((73,T]\). The upper replay
keeps all values below \(1/4\), which verifies the monotonicity domain
for both comparisons. The star family rules out reaching all heads
by uniformly improving this second-moment continuation alone.

Run the [exact verifier](../verify_rh_prime_continuation.py) with Python 3.12 or later from the repository root.

## An unconditional input with the same constants

Both tail calculations use only
\(\lvert\vartheta(x)-x\rvert\le\varepsilon_T x\) for \(x\ge T\).
The RH proof supplies that inequality as above. Published explicit
bounds also give the same two choices of \(\varepsilon_T\) unconditionally.

For \(1423\le x\le10^{19}\), Büthe's Theorem 2, (1.6) and (1.7),
gives \(0<x-\vartheta(x)\le1.95\sqrt{x}\). Hence for
\(T\le x\le10^{19}\),
\[
 \frac{|\vartheta(x)-x|}{x}
 \le \frac{2}{\lfloor\sqrt{T}\rfloor}.
\]
At \(T=200000\), this is at most \(169/(24\cdot447)\).
At \(T=5000000\), it is at most \(256/(24\cdot2236)\).
[Büthe, Theorem 2](https://arxiv.org/abs/1511.02032v2)

For \(x>10^{19}\), use the \(\log X_0=\log X_1=20\) row of
Broadbent–Kadiri–Lumley–Ng–Wilk, Table 14.
Since \(e^{20}<10^{19}\), their lower and upper relative errors
are at most \(9.1639\cdot10^{-5}\) and \(4.2676\cdot10^{-5}\),
respectively. Both are smaller than the two rational
\(\varepsilon_T\) values. The ranges therefore cover every
\(x\ge T\) in each calculation.
[Broadbent–Kadiri–Lumley–Ng–Wilk, Theorem 1 and Table 14](https://arxiv.org/abs/2002.11068)

These estimates can be substituted directly at the partial-summation
step in both continuations. The same rational inequalities then prove
the stated cutoffs, the second-moment threshold, and the applications
to least modulus nine without invoking the full Riemann hypothesis.
