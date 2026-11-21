# Mixed moduli at five and seven when twenty-one is absent

Consider a finite family of residue classes with distinct odd moduli
greater than one, omitting \(3,5,7,21\). If the family covers the integers,
then \(35\) occurs and at least one of \(63,105\) occurs. If \(63\) is
absent, \(15\) occurs and its ternary root has mass greater than \(9/25\)
under the uniform law conditioned to avoid the pure \(3\)-powers.
In particular, when \(27\) is present and the \(15\)- and \(27\)-classes
agree modulo \(3\), modulus \(63\) must occur.

## A common bound with no clipping at five and seven

Use the normalized pure-power reference laws \(\mu_p\) and the clipped
comparison through \(80021\). Take \(\delta_5=\delta_7=0\).
The normalized law starts with the product \(\mu_3\mu_5\mu_7\),
independent of the early mixed classes. Write \(B_5,B_7\) for the
two mixed events under that product, with a class assigned by its
largest prime factor. For any fixed complete terminal query, use the
joint upper envelope \(\Phi_7\) for the conditional suffix cost given
all three coordinates, starting at \(11\). The lower replay from the
preceding result gives \(\Phi_7\ge b_7=81/10000\). After restricting at
\(7\), a conditional upper bound given the coordinates at \(3,5\) is
\[
(1-b_7)\Pr(B_7\mid x_3,x_5)
+\mathbb E[\Phi_7\mid x_3,x_5].
\]
The second term is bounded below by its value when the suffix input
weights with nontrivial \(3\)- or \(5\)-parts are set to zero. The same
lower replay, including the \(7\)-step at zero clipping with the
ordinary reference cap, certifies a lower bound \(b_5=3486/100000\) at
that state. The first term is nonnegative. Restricting at \(5\) then
gives the joint upper bound
\[
\mathbb E\Phi_7 +(1-b_5)\Pr(B_5)+(1-b_7)\Pr(B_7).
\]
This bounds the probability of coverage through \(80021\) plus the
terminal second moment on the survivors, on one normalized law.
[Successive restriction and validation of the floors](82-pure-prime-power-constraints.md)

For \(p=3,5,7\) set \(M_p(0)=1\) and \(M_p(e)=\rho_p p^{-e}\)
when \(e>0\), using the ordinary caps \(\rho_3=6/5\),
\(\rho_5=20/19\), and \(\rho_7=42/41\).
A mixed numerical slot with prime exponents \(e_p\) has mass at most
\(\prod_p M_p(e_p)\). With all slots included, the product caps summed
over the slots coincide with the expected linear charges of the
auxiliary comparison, since \(\Pr(T_p\ge e)=M_p(e)\) for \(e>0\)
and the auxiliary exponents at these three stages are independent.
The upper replay for the full slot sums and the common suffix
expectation gives
\[
U=\frac{108496429444579}{10^{14}}.
\]
The suffix bound is independent of the early mixed slots, so we can
improve \(U\) by subtracting the full slot cap for each known absent
modulus, or the difference between a cap and a better upper bound
on an actual class. Multiply each subtraction by the loss coefficient
\(1-b_p\) at its largest prime. The remaining geometric series are
handled in the replay over every allowed height.

The slot masses for \(21,35,63,105\) are respectively
\(12/205,\ 24/779,\ 4/205,\ 48/3895\). Subtracting the indicated masses
with coefficient \(1-b_7\) proves the first two exclusions in the
table. For the third, let \(\lambda\) be the normalized ternary mass
of the root of the \(15\)-class. If \(15\) is absent, take \(\lambda=0\).
The \(15\)-class has mass at most \((4/19)\lambda\). When
\(\lambda\le9/25\), subtract
\((1-b_5)(4/19)(2/5-9/25)\) as well as the weighted \(21\) and \(63\)
caps.

| Excluded pattern | Joint upper bound |
| --- | ---: |
| \(21,35\) absent | \(77615094537327041/77900000000000000\) |
| \(21,63,105\) absent | \(77535742537327041/77900000000000000\) |
| \(21,63\) absent and \(\lambda\le9/25\) | \(77854834697327041/77900000000000000\) |

All three values are less than \(1-1/2000\). The same unconditional
quadratic continuation leaves an uncovered integer.

If \(15\) and \(27\) agree modulo \(3\), their common root retains Haar
mass at most \(1/3-1/27\). The total remaining ternary mass is at least
\(5/6\), so \(\lambda\le (1/3-1/27)/(5/6)=16/45<9/25\).
Together with the third exclusion this forces \(63\) when \(21\) is
absent. The [exact verifier](../verify_minimum_nine_modulus_fifteen.py)
checks the stated floors and upper bounds.
