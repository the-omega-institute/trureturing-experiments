# Necessary small moduli when the least modulus is nine

A distinct odd cover whose least modulus is \(9\) must include at least
one of \(15,21,33\) and at least one of \(15,21,45\). In particular, if
\(15\) and \(21\) are both absent, then \(33\) and \(45\) must both occur.
The proof below only needs the absence of \(3,5,7\) and one of the two
displayed triples. It permits arbitrary original phases, arbitrary
prime-power heights and any finite prime support.

For context, Dai, Li, Li and Yi prove that the least modulus of a
distinct odd cover belongs to \(\{3,5,7,9\}\) and its least common
multiple is divisible by \(45\). Their bound on the least modulus
makes the case \(9\) the largest remaining minimum. The conclusion
here constrains the numerical moduli in that case.
[Dai–Li–Li–Yi, Theorems 1.1 and 1.2](https://arxiv.org/abs/2611.15110v1)

### Comparison with pure powers removed

For each prime \(p\) normalize Haar on the complement of all actual
pure \(p\)-power classes. Let \(b_p=1\) if \(p\in\{3,5,7\}\) and
\(b_p=0\) otherwise. The least permitted pure exponent is at least
\(b_p+1\). The distinct pure labels therefore delete mass at most
\(\sum_{e\ge b_p+1}p^{-e}=p^{-b_p}/(p-1)\).
Write \(d_p=p-1-p^{-b_p}\). In this normalized reference law every
depth-\(e\) cylinder has mass at most \((p-1)p^{-e}/d_p\).

Through prime \(6997\) use the normalized distortion kernel at
clipping parameter \(\delta_p\) relative to that reference measure.
It preserves the earlier marginal, leaves every pure class outside
its support, and charges the actual forbidden mixed union of
reference fraction \(\alpha\) by \((\alpha-\delta_p)_+/(1-\delta_p)\).
The conditional depth-\(e\) cap after distortion is \(c_pp^{-e}\),
where \(c_p=(p-1)/(d_p(1-\delta_p))\). The chosen parameters give
\(0<c_p\le p\). Primes absent in the original period may be padded to
positive height, with no forbidden classes. As in the preceding
proofs, use all actual heights, including digits used only by later
classes. [Distortion and scalar comparison, Lemmas 2.1, 2.4 and Theorem 2.6](https://arxiv.org/abs/2611.15110v1)

For the old coordinates let \(T_q\) be independent comparison
exponents with \(\Pr(T_q\ge e)=c_q q^{-e}\) for \(e\ge1\).
Write \(D=\prod_{q<p}(1+T_q)\), and let \(\sigma=\min(T_3,2)\)
after processing \(3\), with \(\sigma=0\) initially. Let \(E_p\)
be the set of exponents \(a\) such that \(3^a p\) belongs to the
selected absent triple. The two triples are \(\{15,21,33\}\) and
\(\{15,21,45\}\).

Complete the mixed cofactor inventory with one test phase per
eligible numerical cofactor and current exponent. Distinct full
moduli justify this completion at every exponent. The unit cofactor
is excluded because its class is pure. The additional forbidden
labels in \(E_p\) occur at current exponent one, so their weight is
removed there. Under the convex comparison, the completed mixed
weight at prime \(p\) is
\[
 \frac{D-1-r_p(\sigma)}{d_p},\qquad
 r_p(s)=\frac{p-1}{p}\sum_{a\in E_p}\mathbf1_{\{s\ge a\}}.
\]
Indeed the auxiliary cofactors are divisors of \(\prod_{q<p}q^{T_q}\).
Each nonunit cofactor contributes \(\sum_{e\ge1}(p-1)p^{-e}/d_p=1/d_p\).
Each missing cofactor \(3^a\) at current exponent one removes
\((p-1)/(pd_p)\) exactly when \(T_3\ge a\). This comparison uses
nonnegative weights on the permitted labels before evaluating the
positive part, so the subtraction represents absent slots in the
same comparison and the same original family. Put
\(U_p(s)=1+\delta_p d_p+r_p(s)\). The actual stage loss is at most
\[
 \frac{\mathbb E(D-U_p(\sigma))_+}{d_p(1-\delta_p)} .
\]
For every complete divisor query, the same comparison applied to
the square bounds its second moment at the end of the prefix by
\[
 K=\prod_{3\le p\le6997,\ p\ {\rm prime}}
       \left(1+c_p\frac{3p-1}{(p-1)^2}\right).
\]
All these estimates hold on the same distortion law. Restrict it
once to the event that survives every stage through \(6997\).
The restricted measure has mass at least one minus the sum of
stage losses, and its second moment is still bounded by \(K\).

### Exact prefix calculation and infinite continuation

Use the following \(\delta_p\) in hundredths. The first row lists
inclusive upper endpoints, with each interval starting just above
the previous endpoint.

| Upper endpoint | 3 | 5 | 7 | 11 | 13 | 17 | 19 | 23 | 29 | 31 | 37 |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| \(100\delta_p\) | 0 | 0 | 0 | 12 | 18 | 20 | 18 | 24 | 26 | 24 | 20 |

| Upper endpoint | 41 | 43 | 47 | 53 | 101 | 211 | 419 | 839 | 1679 | 3359 | 6997 |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| \(100\delta_p\) | 28 | 26 | 26 | 28 | 28 | 30 | 30 | 30 | 30 | 28 | 26 |

Keep all states \((j,s)\) with \(1\le j\le1850\) and \(0\le s\le2\),
using \(S=10^{24}\) and upper integer values for
\(S\Pr(D=j,\sigma=s)\). Start with mass \(S\) at \((1,0)\).
At a prime \(p\), each multiplier \(i=1+T_p\) maps \(j\) to \(ij\).
At \(3\) set \(s\) to \(\min(i-1,2)\). At other primes retain \(s\).
Multiply the stored integers by the ceilings of
\(S\Pr(1+T_p=i)\), sum by destination, divide by \(S\), and round each
state upward. The multiplier masses at \(i=1\) and \(i\ge2\) are
\(1-c_p/p\) and \(c_p(p-1)/p^i\), respectively. Larger states cannot
return to the retained range because every multiplier is at least one.

Independently carry upper bounds for the full mean of \(D\) and its
second moment, starting each at \(S\). The update factors are
\(1+c_p/(p-1)\) and \(1+c_p(3p-1)/(p-1)^2\), with rounding upward.
Before each update, use the exact loss numerator
\[
 \mathbb E(D-U_p(\sigma))_+
 =\mathbb E D-\mathbb E U_p(\sigma)
   +\sum_{j,s}(U_p(s)-j)_+\Pr(D=j,\sigma=s).
\]
Every positive term in the sum is retained. The expectation of
\(U_p(\sigma)\) uses the exact probabilities
\(\Pr(T_3\ge a)=6/(5\cdot3^a)\) for the only nonzero \(E_p\), because
\(\delta_3=0\) and \(c_3=6/5\). Thus all terms replaced by upper bounds
have nonnegative coefficients. Round each stage loss up at scale \(S\)
before summing.

Apply the [second-moment continuation](76-rh-prime-weighted-continuation.md#the-second-moment-threshold-at-73) to the restricted
measure at \(6997\), using RH for the infinite remainder. Its same
upper replay, starting at \(5000000\) and stopping before processing
\(6997\) in decreasing order, gives
\[
 v=\frac{27144070771294538066}{S}.
\]
The finite prefix passes if its remaining mass is greater than \(Kv\).
The following table gives the exact scaled integers for the prefix
bounds and the resulting strict mass margin.

| Absent moduli in addition to \(3,5,7\) | \(S\) times the sum of loss bounds | \(S\) times the moment bound | Remaining mass minus \(Kv\), greater than |
| --- | ---: | ---: | ---: |
| \(15,21,33\) | 938773600338587870460242 | 1908984849362909824839827392 | \(9/1000\) |
| \(15,21,45\) | 935988366744625136574530 | 1908984849362909824839827392 | \(12/1000\) |

Every later class is processed by the standard continuation, including
pure powers of larger primes. The final positive live mass gives an
integer outside the original classes by CRT. This rules out the two
missing-label patterns while leaving open the case with least modulus
\(9\) and either \(15\) or \(21\), and the case with both \(33\) and
\(45\) and neither \(15\) nor \(21\).

Run the [exact verifier](../verify_rh_prime_continuation.py) with Python 3.12 or later from the repository root.
