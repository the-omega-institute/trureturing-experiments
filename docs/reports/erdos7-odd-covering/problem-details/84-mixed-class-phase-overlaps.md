# Relative phases of fifteen, thirty-five and sixty-three

Consider a finite covering by distinct odd moduli greater than one,
with \(3,5,7,21\) all absent. The classes of moduli \(15\) and \(35\)
occur. If they agree modulo \(5\), then \(63\) occurs. In this case,
if \(63\) and \(15\) also agree modulo \(3\), then \(105\) occurs with
a class disjoint from each of the \(15\)- and \(35\)-classes.
These statements impose no prime-power height or support bound.
[Presence of the early moduli](83-early-mixed-modulus-constraints.md)

## Charging the late classes only on the early survivors

Use the product reference law at \(3,5,7\), zero clipping at \(5,7\),
and the comparison at \(11\) and above from the preceding proof.
Write \(C_m\) for the original class of modulus \(m\), and let \(B_5\)
and \(B_7\) be the mixed unions at their respective primes. Write
\(b_5,b_7\) for the previously certified floors and \(w_p=1-b_p\).
For a selected set \(\mathcal{E}\) of actual classes assigned to \(7\),
bound the conditional probability of \(B_7\) by the probability of its
remaining union \(B_7'\) plus the sum of the probabilities of \(C\in\mathcal{E}\).
As in the preceding proof, conditional on \(x_3,x_5\) this bounds the
suffix by
\[
J(x_3,x_5)+w_7\sum_{C\in\mathcal{E}}\Pr(C\mid x_3,x_5),
\qquad
J=\mathbb E[\Phi_7\mid x_3,x_5]+w_7\Pr(B_7'\mid x_3,x_5).
\]
The lower replay at zero clipping gives \(\mathbb E[\Phi_7\mid
x_3,x_5]\ge b_5\) through the state with all nontrivial processed
weights set to zero, so \(J\ge b_5\).
Restricting at \(5\) applies the \(b_5\) credit to \(J\) and
integrates the selected terms on \(B_5^{\mathrm c}\).
Since \(C_{15}\subseteq B_5\), the resulting upper bound is
\[
\mathbb E\Phi_7+w_5\Pr(B_5)+w_7\Pr(B_7')
+w_7\sum_{C\in\mathcal{E}}\Pr(C\setminus C_{15}).
\]
It uses the same suffix and product law as the full-slot bound \(U\).
In \(U\), replace the affected terms by a joint upper bound on
\(w_5\Pr(C_{15})+w_7\sum_{C\in\mathcal{E}}\Pr(C\setminus C_{15})\).
All earlier coefficients and the continuation are as in the preceding
proof. [Full-slot upper bound and floors](83-early-mixed-modulus-constraints.md)

## Exact corrections for shared residues

Let \(\lambda\) be the ternary mass of the root of \(C_{15}\) and let
\(s\) be the quinary mass of its root. If \(C_{35}\) has the same
quinary root, its contribution outside \(C_{15}\) is at most
\(s(6/41)(1-\lambda)\). The joint charge of \(15,35\) is at most
\[
w_5 s\lambda + w_7 s(6/41)(1-\lambda).
\]
This expression increases with \(s\) and \(\lambda\), since
\(w_5>w_7(6/41)\). Take \(s=4/19\) and \(\lambda=2/5\).
The resulting correction to the full-slot bound is a decrease of
\(w_7(48/3895)\). If \(63\) is absent, the weighted \(21,63\) caps
are also omitted, giving the bound
\[
U-w_7\left(\frac{12}{205}+\frac4{205}+\frac{48}{3895}\right)
=\frac{77535742537327041}{77900000000000000}.
\]
This also excludes \(63\) being contained in a \(9\)-class, since
deleting a contained class preserves coverage.

When \(C_{63}\) has the same ternary root as \(C_{15}\), let \(\tau\)
be the mass of the finer ternary cylinder of \(C_{63}\). Its mass
outside \(C_{15}\) is at most \(\tau(6/41)(1-s)\).
The upper bound for the three terms together is then
\[
w_5 s\lambda + w_7(6/41)\bigl(s(1-\lambda)+\tau(1-s)\bigr).
\]
It increases with \(\lambda,\tau\) and with \(s\) on the whole domain,
because \(\lambda+\tau\le2/5+2/15<1\).
Use \(s=4/19\), \(\lambda=2/5\), and \(\tau=2/15\).
The correction to \(U\) is a decrease of
\(w_7(48/3895+(4/205)(4/19))\). If \(105\) is absent, the additional
savings from omitting \(21,105\) give
\[
U-w_7\left(\frac{12}{205}+\frac{48}{3895}+
                 \frac{48}{3895}+\frac4{205}\frac4{19}\right)
=\frac{77773798537327041}{77900000000000000}.
\]
Since both \(15\) and \(35\) divide \(105\), a \(105\)-class that
intersects either class is contained in it and could be deleted.

Both rational bounds are below \(1-1/1000\), as checked by the
[exact verifier](../verify_minimum_nine_modulus_fifteen.py).
The quadratic continuation therefore leaves an uncovered integer
in each excluded case.
