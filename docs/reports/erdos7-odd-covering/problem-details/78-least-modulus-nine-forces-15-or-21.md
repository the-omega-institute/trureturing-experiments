# Least modulus nine forces 15 or 21

The extended joint comparison proves that [least modulus nine forces 15](79-least-modulus-nine-forces-15.md).

Every finite distinct odd cover with least modulus \(9\) contains a
class of modulus \(15\) or \(21\). More generally, a distinct odd family
with no moduli \(3,5,7,15,21\) leaves an integer uncovered. The residues,
prime-power heights and finite prime support are arbitrary.

The proof uses the pure-power reference measures and prime schedule
from the preceding report, and the same RH continuation beyond \(6997\).
Here the conditional density bound depends on the actual forbidden fibre
at every prime. Its contribution is propagated through the later loss
estimates and the second moment used at the end of the prefix.
[Pure-power reference and schedule](77-necessary-small-moduli-at-minimum-nine.md) ·
[RH continuation](76-rh-prime-weighted-continuation.md#the-second-moment-threshold-at-73)

## A joint bound for the prefix losses and terminal moment

Use one normalized distortion law on the prime coordinates through
\(6997\). At each prime \(p\), let \(\alpha\) be the mixed forbidden
fraction under the reference law conditioned away from the actual pure
powers. Write \(\rho=(p-1)/d_p\), \(\delta=\delta_p\), and
\(g=(1-\delta)^{-1}\). The conditional cap for a depth \(e>0\) cylinder
is \(\rho p^{-e}/(1-\min(\alpha,\delta))\). It is at most
\(\bar c p^{-e}\), where \(\bar c=\rho g\).

Fix an arbitrary complete divisor query on the entire prefix, with
load \(L\). The quantity to bound above is the sum of the assigned
mixed-event probabilities through \(6997\), plus \(v\mathbb E L^2\)
on the same normalized law. Here
\(v=27144070771294538066/10^{24}\) is the second-moment continuation
threshold for all later primes. If the bound is below one for every
query, restriction to the common survivors has mass greater than \(v\)
times the maximum complete-query second moment on that restriction.

We use the joint convex comparison of Dai, Li, Li and Yi with the
nonuniform reference law. The following normalization gives the needed
form of their Lemma 3.2. For an upper bound \(a\) on \(\alpha\), put
\(c(a)=\rho/(1-\min(a,\delta))\). For a directionally convex increasing
continuation \(\Phi\), write \(V_c\) for its nested-prefix comparison
with exponent tails \(c p^{-e}\). The dependence on \(c\) is affine.
Dai–Li–Li–Yi's upper envelope becomes
\[
 g(a-\delta)_+ +
 \max\{V_{c(a)},\ V_{\bar c}-(\bar c-c(a))/\bar c\}.
\]
The first branch bounds the conditional continuation using the actual
cap. For the directional convexity check use their two-branch argument
with coefficient bound \(1/\bar c\). At \(a=\delta\), the left derivative
in the fibre fraction is at most
\(\rho g^2/\bar c=g\), the right derivative of the clipping loss.
The coefficient of the cap in \(V_c\) is nonnegative and increasing
in the other weights, as in their equation (31). This gives the same
increasing differences argument for the joint envelope. The partition
into present cofactor slots and omitted numerical labels uses
nonnegative weights on the present slots, so the linear substitutions
in their backward comparison still have nonnegative coefficients.
[Directional convex comparison, Lemmas 3.1 and 3.2](https://arxiv.org/html/2611.15110v1#S3)

Apply the comparison backwards from the terminal test \(vL^2\).
At the terminal step the square is increasing and directionally
convex in all the nonnegative contributions to \(L\).
Each comparison uses the same auxiliary exponent for all later tests
and the terminal query. A complete fixed query contributes one term
for each divisor. In the auxiliary layout its load is the product
of the surviving depth counts \(D\). At prime \(p\) the completed
mixed weight, before the positive part, is
\((D-1-r)/d_p\). The unit is absent because its label is pure.
At \(p=5,7\), the missing \(15,21\) removes
\(r=(p-1)/p\) exactly when the auxiliary ternary exponent is positive.
At every other prime \(r=0\). These are precisely the weights of the
admitted divisor slots under the same joint comparison. First use
finite prefix windows at each auxiliary prime, then enlarge to all
heights by the monotone comparison of the two-branch operators.
Every suffix expectation is bounded by a quadratic in \(D\), as
verified below using geometric tails. This supplies a finite upper
bound throughout the enlargement.

## Scalar replay of the joint bound

For \(p\ge11\) the state is \(D\). At \(5\) and \(7\) also retain
whether \(T_3>0\). Let \(\Phi(D)\) denote the continuation after the
current prime. For \(N=1+T_p\) with cap \(c\), write
\(V_c(D)=\mathbb E_c\Phi(DN)\). Then use the displayed joint envelope
with \(a=\max(0,(D-1-r)/d_p)\). Start beyond \(6997\) with
\(\Phi(D)=vD^2\). Above the state cutoff \(W=1850\), all clipping
thresholds \(1+r+\delta d_p\) are below \(D\), so the state value is
exactly a quadratic \(AD^2+BD+C\). Given the coefficients for the
continuation, the new high-state coefficients are
\[
 A'=A\left(1+\bar c(3p-1)/(p-1)^2\right),\qquad
 B'=B\left(1+\bar c/(p-1)\right)+g/d_p,\qquad
 C'=C-g(1+r+\delta d_p)/d_p.
\]
For \(1\le D\le W\), put \(m=\lfloor W/D\rfloor\).
Use the point masses \(1-c/p\) for \(N=1\) and \(c(p-1)/p^i\)
for \(2\le i\le m\). On \(N>m\), of probability \(c/p^m\),
write \(N=m+G\), where \(G\) is geometric on positive integers.
Its first two moments are \(p/(p-1)\) and \(p(p+1)/(p-1)^2\).
Use these full tail moments with the high-state polynomial.
At \(3\), use \(D=1\) with cap \(6/5\) and zero mixed loss.
The \(N=1\) branch has \(T_3=0\). The rest has \(T_3>0\).

The replay rounds all table values and all three polynomial
coefficients upward on the grid \(10^{-24}\). For the strictly
positive-depth expectation per unit cap, round its point probabilities
and its remaining tail probability upward before evaluating the
expectation. This is valid because the computed high-state polynomial
bounds a nonnegative continuation, and the coefficient of \(N=1\)
in the full expectation is nonnegative. Both branches of the
joint envelope are increasing in the stored continuation values.
The signed constant \(C\) is kept with its sign in each tail expression.

The upper bound after the 896 stages from \(6997\) through \(11\) and
the three smaller-prime steps is
\[
 \frac{986888733243508975907980}{10^{24}} < \frac{987}{1000}.
\]
Thus the prefix survivor mass exceeds its required second-moment
reserve by more than \(13/1000\). All later primes are handled on
that one survivor measure by the standard continuation. A positive
final mass yields an uncovered integer.

Run the [exact verifier](../verify_joint_clipping_minimum_nine.py)
with Python 3.12 or later from the repository root.

## A phase constraint from the modulus-nine class

Suppose a distinct odd cover has least modulus \(9\). A modulus in this
section is **inactive** if it is absent or its whole class is contained
in the class of modulus \(9\). Then both of the following constraints
hold.

- If \(15\) is absent and the \(21\)-class agrees with the \(9\)-class
  modulo \(3\), then \(45\) is present and its class is active.
- If \(21\) is absent and the \(15\)-class agrees with the \(9\)-class
  modulo \(3\), then at least one of \(45\) and \(63\) is present with
  an active class.

Delete any classes contained in the \(9\)-class before running the same
proof. This preserves the union and the pure-power reference laws at
\(3,5,7\). Let \(\lambda\) be the normalized Haar law on the surviving
ternary coordinate after removing every pure \(3\)-power class.
For a present class of modulus \(j\in\{15,21\}\), let \(\lambda_j\)
be the mass of its ternary residue modulo \(3\). Set \(\lambda_j=0\)
when the class is absent. On the root specified by the \(9\)-class, the
pure survivors have Haar mass at most \(1/3-1/9=2/9\).
The total surviving Haar mass is at least \(5/6\), so this
root has \(\lambda\)-mass at most \(4/15\).

The clipping parameters at \(3,5,7\) are zero. Under the original
normalized distortion law their marginal is consequently the
product of the three pure-power reference laws. At \(5\) and \(7\)
the charged probability is bounded above by the sum of the individual
mixed class probabilities under this same product. The two slots
\(15\) and \(21\) were excluded in the preceding calculation.
If present, their contributions here are at most
\((4/19)\lambda_{15}\) and \((6/41)\lambda_{21}\), respectively,
using the cap at the other prime.

The completed ledger at \(5\) included the \(45\)-slot with
coefficient \((4/19)(6/5)3^{-2}=8/285\). At \(7\) it included the
\(63\)-slot with coefficient \((6/41)(6/5)3^{-2}=4/205\).
If either is inactive, remove its coefficient from the corresponding
completed upper bound. Each subtraction removes a term from the
completed sum of cylinder caps. The bound for the joint suffix starting
at \(11\), including the terminal query, uses the same three-coordinate
product and the same caps as above. In particular, it is independent
of which of these four early slots occurred. The early loss sum is
linear in the admitted slots since the clipping parameters there are zero.

Let \(C_0=986888733243508975907980/10^{24}\) be the outward-rounded
bound above with both \(15\) and \(21\) absent and both \(45\) and
\(63\) completed. The exact auxiliary upper expression computed before
rounding is the sum of the two completed early losses and the joint
suffix evaluated on that product comparison law. Each changed slot
adds or removes the indicated expected contribution before rounding.
The unchanged suffix permits the same \(C_0\) to serve in the
resulting bound
\[
 C_0+\frac4{19}\lambda_{15}+\frac6{41}\lambda_{21}
  -\frac8{285}\mathbf 1_{\{45\ {\rm inactive}\}}
  -\frac4{205}\mathbf 1_{\{63\ {\rm inactive}\}}.
\]
If this expression is below one, restriction to survivors followed
by the same RH continuation leaves an uncovered integer.

Under the first set of excluded conditions, this bound is at most
\(C_0+(6/41)(4/15)-8/285<1-2/1000\).
Under the second it is at most
\(C_0+(4/19)(4/15)-8/285-4/205<1-4/1000\).
The same verifier checks both strict margins.
