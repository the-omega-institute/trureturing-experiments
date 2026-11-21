# Least modulus nine forces modulus fifteen

Every finite covering of the integers by residue classes with distinct odd
moduli and least modulus \(9\) contains a class of modulus \(15\).
More generally, every finite family with distinct odd moduli greater than
one and with \(3,5,7,15\) all absent leaves an integer uncovered. The
argument allows \(21\) and arbitrary original residues, prime-power
exponents and finite prime support.

## Joint comparison through prime \(80021\)

Use the normalized pure-power reference measures and the joint
comparison for mixed losses and the terminal second moment from the
preceding report. The conditional reference cap at each prime is
\(\rho p^{-e}\) at depth \(e>0\), where
\[
 d_p=p-1-p^{-b_p},\qquad
 \rho=(p-1)/d_p,\qquad
 b_p=\mathbf 1_{\{3,5,7\}}(p).
\]
These caps apply because the moduli \(3,5,7\) are absent. Process primes
through \(80021\) with the schedule in the verifier, which uses positive
clipping already at \(7\). As before, pad any missing prime to a
coordinate with no assigned classes. All terminal queries are complete
on the actual heights, including heights from later classes.

At a given step let \(D\) be the product of the earlier auxiliary depth
counts. In the mixed-fibre comparison put
\[
 a=\frac{(D-1-r)_+}{d_p}.
\]
The unit cofactor belongs to the pure-power classes already removed.
At \(p=5\) the missing modulus \(15\) removes the cofactor \(3\) at
depth one. Thus \(r=4/5\) at \(p=5\) when the auxiliary ternary depth
is positive, and \(r=0\) in every other state and stage. All other
mixed slots enter the complete upper comparison, including \(21\)
at \(p=7\). This assigns a nonnegative weight to each admitted slot
before applying the joint convex comparison.
[Joint clipping bound and same-law argument](78-least-modulus-nine-forces-15-or-21.md)

For a fixed complete terminal query with load \(L\), charge all
assigned mixed losses and the term \(v\mathbb E L^2\) on the same
normalized law. Here \(v\) is the upper continuation threshold beyond
\(80021\) obtained by the exact quadratic replay from the prime tail.
The replay captures
\[
 v=\frac{1697239918308703913}{10^{24}}
\]
before processing \(80021\) in decreasing order. The analytic input to
this replay can use the unconditional Chebyshev bounds.
[Quadratic replay and unconditional prime estimates](76-rh-prime-weighted-continuation.md)

Starting with the terminal value \(vD^2\), apply the joint envelope
at each prime in reverse order. With \(g=(1-\delta)^{-1}\),
\(\bar c=\rho g\), \(c(a)=\rho/(1-\min(a,\delta))\), and
\(V_c(D)=\mathbb E_c\Phi(DN)\) for \(N=1+T_p\), the update is
\[
 g(a-\delta)_+
 + \max\{V_{c(a)}(D),\,V_{\bar c}(D)-(\bar c-c(a))/\bar c\}.
\]
The auxiliary exponent is shared by the later mixed-loss tests
and the terminal query, and
\(\Pr_c(T_p\ge e)=cp^{-e}\) for \(e\ge1\).
The \(5\)-step uses two separate updates for zero and positive ternary
depth. At \(3\) use \(D=1\), no mixed loss and the cap \(6/5\).
This puts mass \(3/5\) on zero ternary depth.

## Exact finite replay with geometric tails

The table stores an upper value for every \(1\le D\le W=34500\).
Every clipping threshold \(1+r+\delta d_p\) is at most \(W\).
On \(D>W\), the same high-state quadratic \(AD^2+BD+C\) from
the preceding report is an upper value. The verifier rounds all its
coefficients and all table values upward in units of \(10^{-14}\).
It starts the quadratic coefficient at
\(169723992/10^{14}\), the upward rounding of \(v\).

The function \(AD^2+BD\) with those same nonnegative coefficients
also bounds the continuation at every positive integer \(D\).
For the induction, dominate the mixed charge by \(gD/d_p\), and the
future expectation by \(V_{\bar c}\) applied to this increasing
polynomial. The resulting coefficients are exactly the updates of
\(A,B\) in the high-state quadratic. This permits a short finite sum
for each expectation without discarding any mass.

In the expectation per unit positive-depth cap, sum the point terms
for \(2\le N\le m\), with
\[
 m=\min\bigl(\lfloor W/D\rfloor,\ \min\{j\ge1:p^j\ge10^{14}\}\bigr).
\]
Round those point probabilities upward on the same grid. On \(N>m\)
use its full geometric mass per unit cap \(p^{-m}\) and write
\(N=m+G\), with \(G\) positive geometric of mean \(p/(p-1)\) and
variance \(p/(p-1)^2\). If \(m=\lfloor W/D\rfloor\), all remaining
states are above \(W\), so evaluate \(AD^2+BD+C\) at the random
argument \(DN\) using these two moments. Otherwise evaluate the
global upper bound \(AD^2+BD\) at \(DN\). The verifier takes the
ceiling of the rational tail expression in both cases. Each
\(V_c(D)\) combines this positive-depth expectation per unit cap
with the value for \(N=1\), whose coefficient \(1-c/p\) is nonnegative.
The subtracted difference of caps in the second branch of the envelope
is rounded downward.

The replay covers 7834 stages at primes from \(11\) through \(80021\)
and the three smaller-prime steps. Its joint upper bound is
\[
 \frac{99718016274333}{10^{14}} < 1-\frac1{500}.
\]
The estimate holds for every fixed complete query. Restrict the normalized
law to its actual survivors through \(80021\). Its mass then exceeds
\(v\) times its maximum complete-query second moment. The quadratic
continuation gives positive survivor mass through all remaining primes.
An integer in the surviving residue set misses every class.

Run the [verifier](../verify_minimum_nine_modulus_fifteen.py) from the
repository root with Python 3.12 or later and a GCC or Clang C++17
compiler with signed 128-bit integers and undefined-behavior sanitization.
It computes the continuation, compiles the integer prefix replay in
a temporary directory, and verifies the strict reserve.
