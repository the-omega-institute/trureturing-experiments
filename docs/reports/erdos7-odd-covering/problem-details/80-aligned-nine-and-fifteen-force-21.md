# Residues of the moduli nine and fifteen when twenty-one is absent

Suppose a finite covering by distinct odd moduli has least modulus \(9\).
It contains modulus \(15\). If modulus \(21\) is absent, the classes of
moduli \(9\) and \(15\) have different residues modulo \(3\), and modulus
\(45\) occurs with a class outside the class of modulus \(9\). The proof
allows arbitrary residues and prime-power exponents.
[Absence of modulus \(21\) also forces modulus \(33\)](81-twenty-one-or-thirty-three.md)
[Presence of modulus \(15\)](79-least-modulus-nine-forces-15.md)

## Restricting at prime five

Use the joint comparison and quadratic continuation of the preceding
result, changing the clipping parameter at \(7\) to \(3/100\). Keep
\(\delta_5=0\). Suppose \(21\) is absent. At \(7\) omit the depth-one
slot with ternary cofactor \(3\), reducing the completed mixed weight by
\((6/7)/d_7\) when the auxiliary ternary depth is positive. At \(5\)
handle the class of modulus \(15\) separately and omit its slot in the
remaining completed mixed weight. The normalized distortion at \(5\)
equals the reference law because \(\delta_5=0\), regardless of the
forbidden mixed set.

Fix a complete terminal query \(L\). On the normalized law write \(B_p\)
for the mixed forbidden set assigned to prime \(p\), and set
\[
S=\sum_{p\ge7}\mathbf 1_{B_p}+vL^2 ,
\]
where the sum ends at \(80021\) and \(v\) is the continuation threshold.
The pointwise upper bound
\[
\mathbf 1_{\mathrm{prefix\ covered}}+
 vL^2\mathbf 1_{\mathrm{prefix\ survivor}}
 \le \mathbf 1_{B_5}+\mathbf 1_{B_5^{\mathrm c}}S
\]
counts the \(5\)-stage loss before the suffix. The pure powers have
already been removed by the reference laws. Let \(\Phi(z)\) be the joint
upper bound for the conditional expectation of \(S\) given the ternary
and quinary coordinates, where \(z\) is the vector of surviving cofactor
weights. The event \(B_5\) is determined by these two coordinates. If
\(\Phi(z)\ge b\) everywhere, then the expectation on the right is at most
\[
\mathbb E\Phi(z)+(1-b)\Pr(B_5).
\]
Indeed, subtracting the expectation of \(\Phi(z)\) on \(B_5\) saves at
least \(b\Pr(B_5)\).

To bound \(\Phi(z)\) from below for this subtraction, complete the
suffix inventories with arbitrary phase tests in the admissible slots.
The bound is increasing in their nonnegative weights by the joint
comparison. The coordinates of \(z\) corresponding to suffix weights
with trivial \(3\)- and \(5\)-parts retain their full values at every
point. Setting every other coordinate to zero therefore gives a
pointwise lower bound on \(\Phi\). Its scalar state immediately after
\(5\) is \(D=1\) with ternary depth zero. In this state the omission of
\(21\) has no effect at \(7\). One may extend the finite periods by
adding digits, leaving the original family unchanged. Take the window
wide enough to contain every retained exponent and completed slot used
in the floor computation. The upper comparison allows all exponents as
in the preceding report.

For a floor on this value, keep only the primes from \(7\) to \(211\),
start with zero after \(211\), and use the first branch of the joint
envelope at each step. Retain only \(1\le D\le400\) and discard the
remaining geometric probabilities in the positive-depth expectation.
The first branch is the clipping charge plus the expectation with cap
\(c(a)\). Both terms are nonnegative, so downward rounding of each
term gives a lower bound. Monotonicity in the continuation propagates
the bound to the initial state \(D=1\). The downward replay yields
\[
b=\frac{3}{100} \ \le\ \frac{3081633670362}{10^{14}} .
\]

## Exact upper bound

Without either slot \(15\) or \(21\), the unweighted joint comparison
with the new \(7\)-stage clipping gives the upper value
\[
U_0=\frac{94667973832900}{10^{14}} .
\]
At \(5\) the expectation of the completed early charge for the other
mixed slots, in this scalar comparison, is exactly
\[
\frac{\mathbb E T_3-(4/5)\Pr(T_3\ge1)}{d_5}
=\frac{3/5-(4/5)(2/5)}{19/5}=\frac7{95}.
\]
The expectation is linear in the coefficient of this charge because
\(\delta_5=0\). In particular, the joint directional comparison for
\(\mathbb E\Phi(z)\) plus the charge weighted by \(1-b\) is bounded by
\(U_0-b(7/95)\). The \(5\)-coordinate law, hence the suffix expectation
and its comparison, is independent of the phase of \(15\).

If the \(15\)-class agrees with the \(9\)-class modulo \(3\), their
common root has normalized ternary mass at most
\((1/3-1/9)/(5/6)=4/15\). The depth-one cylinder on the
\(5\)-coordinate has reference mass at most \(4/19\). The contribution
of the \(15\)-class on the product of the two reference laws is
therefore at most \(16/285\). The union bound at \(5\) and the same
charge coefficient \(1-b\) give the final bound
\[
U_0-\frac3{100}\frac7{95}+\frac{97}{100}\frac{16}{285}
<\frac{999}{1000}.
\]
This uniformly bounds the cover charge plus the terminal moment
restricted to survivors, for every complete query. The quadratic
continuation then gives an uncovered integer.

The same bound also applies when \(45\) is absent or its class is
contained in the class of modulus \(9\). In the latter case delete it
without changing the union. Omitting the slot with ternary cofactor \(9\)
at depth one at \(5\) saves
\[
\frac4{19}\Pr(T_3\ge2)=\frac8{285}
\]
from the unweighted early comparison. The suffix expectation remains
the same because \(\delta_5=0\). With arbitrary phase, the contribution
of the \(15\)-class is bounded by
\((4/19)(2/5)=24/285\). Their net contribution is at most \(16/285\),
with the same common coefficient \(1-b\).

The [exact verifier](../verify_minimum_nine_modulus_fifteen.py) checks
both rounded replays and the strict rational inequality. The same
unconditional prime estimates as in the preceding result certify the
continuation.
