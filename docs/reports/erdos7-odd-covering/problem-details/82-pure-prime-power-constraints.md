# Pure prime powers in covers with no modulus below nine

Let a finite family of residue classes have distinct odd moduli greater
than one, with \(3,5,7\) absent. Each of the following conditions suffices
for an uncovered integer.

- Modulus \(21\) is absent and the union of the pure \(3\)-power classes
  has uniform measure at most \(5/36\).
- Moduli \(21\) and \(25\) are absent.
- Moduli \(25,27,49\) are all absent.

There are no further restrictions on phases, exponents or the finite set of
prime factors. In the first condition the uniform measure is taken modulo
the largest pure \(3\)-power modulus, with zero measure for an empty
family of pure powers.

In particular, a cover of least modulus \(9\) without \(21\) must contain
\(25\) and a \(27\)-class disjoint from its \(9\)-class. If \(27\) is absent,
the pure ternary union has measure at most
\[
3^{-2}+\sum_{e\ge4}3^{-e}=\frac7{54}<\frac5{36}.
\]
A \(27\)-class contained in the \(9\)-class can be deleted without changing
the union. The same deletion argument requires the \(45\)-class from the
earlier phase constraint to be disjoint from the \(15\)-class as well.
[Earlier phase constraint](80-aligned-nine-and-fifteen-force-21.md)

## Reference measures

Remove the pure power classes at each prime from its uniform reference
measure and normalize. If the union at \(p\) has uniform measure at
most \(h_p\), a depth-\(e\) cylinder has reference mass at most
\(\rho_p p^{-e}\), where \(\rho_p=(1-h_p)^{-1}\). Write
\(d_p=(p-1)/\rho_p\). The complete inventory at a prime has mixed
weight \((D-1)/d_p\). The unit cofactor is assigned to pure powers.

At \(3\) use auxiliary cap \(36/31\) when \(h_3\le5/36\), cap \(54/47\)
when \(27\) is absent, and cap \(6/5\) otherwise. With the square of
\(p=5\) or \(p=7\) absent, pure powers at that prime have exponents
at least three. Their union has measure at most \(p^{-2}/(p-1)\),
so use \(d_p=p-1-p^{-2}\) for that prime. All other caps are those
of the [joint comparison](79-least-modulus-nine-forces-15.md).

When \(21\) is absent, clip at \(7\) with \(\delta_7=3/100\) and omit
the depth-one slot with cofactor \(3\) from its mixed inventory.
At \(5\) use \(\delta_5=0\) and charge \(97/100\) of the mixed loss.
The previous lower calculation through \(211\) supplies the
\(3/100\) credit on the pointwise suffix bound. The change of the
ternary or quinary cap leaves the floor calculation at its initial
state \(D=1\) independent of both caps. This gives the first two bounds
in the table below. [Prime-five credit](80-aligned-nine-and-fifteen-force-21.md)

## Coupling the losses at five and seven

For the third condition use \(\delta_7=17/100\) and allow all mixed slots.
In this case charge \(1-b_p\) of the mixed loss at each of \(p=7,5\),
with
\[
b_7=\frac{81}{10000},\qquad b_5=\frac{3486}{100000}.
\]
At a stage on the normalized law, let \(B_p\) be the event removed when
restricting to surviving points. Conditional on the earlier coordinates,
its probability is \(g(a)=(a-\delta_p)_+/(1-\delta_p)\), where \(a\) is
the actual forbidden reference fraction. Suppose a pointwise upper
bound \(\Phi\) for the suffix cost is everywhere at least \(b_p\).
The suffix includes the terminal second moment only on surviving
points. Replacing the suffix by \(\Phi\) on the complement of \(B_p\)
bounds the cost by
\[
\Pr(B_p)+\mathbb E[\boldsymbol{1}_{B_p^{\mathrm c}}\Phi]
\le (1-b_p)\Pr(B_p)+\mathbb E\Phi .
\]
This inequality holds conditionally on each earlier point. Apply it at
\(7\) before applying it at \(5\). All costs refer to one normalized law.

For a stage with weight \(w=1-b_p>0\), apply the joint comparison with
suffix \(\Phi/w\) and multiply by \(w\). The resulting envelope is
\[
w g(a)+\max\!\left\{V_{c(a)},\
V_{\bar c}-w\frac{\bar c-c(a)}{\bar c}\right\},
\]
where \(a\) can now be the completed upper bound on the forbidden
fraction, and \(V_c\) denotes the expectation of \(\Phi\) with cap
\(c\). The same positive cap kernels and increasing directional
convexity apply. This also gives the changes to the high-state quadratic
and the signed constant by multiplying the current clipping charge
by \(w\). [Joint comparison lemma](https://arxiv.org/html/2611.15110v1#S3.Thmtheorem2)

To verify each floor, complete the present inventories and set all
suffix input weights with nontrivial processed prime parts to zero.
Monotonicity makes this a lower bound at every actual point, with
initial state \(D=1\). At \(7\) this floor uses only primes above \(7\)
and the terminal query \(vD^2\), with \(v=169723992/10^{14}\) rounded
up from the valid continuation threshold. To check the floor at \(5\),
include the already weighted \(7\)-step. In both cases retain states
\(1\le D\le1200\), drop the other transitions, round retained terms
down, and round the subtracted cap difference up. The maximum of the
two resulting branches is still a lower bound. Finite prime-power
windows can include all the retained auxiliary states by adding digits.
The upper replay still uses every exponent, with its geometric tails
as before.

## Rational bounds

The upper and lower replays use units of \(10^{-14}\). Their upper
values, after the credits, are as follows.

| Condition | Upper numerator |
| --- | ---: |
| \(21\) absent and the ternary union has measure at most \(5/36\) | \(99906583004156\) |
| \(21,25\) absent | \(99995591954516\) |
| \(25,27,49\) absent | \(99995695893446\) |

The first bound leaves more than \(1/2000\) of reserve and the other
two leave more than \(1/50000\). For each complete query, these bounds
control the union of covered points plus the terminal moment on its
complement. The unconditional quadratic continuation therefore gives
an uncovered integer. The [exact verifier](../verify_minimum_nine_modulus_fifteen.py)
checks every floor and upper bound.
