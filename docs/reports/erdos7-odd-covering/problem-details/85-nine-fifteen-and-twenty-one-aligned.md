# Pure powers when nine, fifteen and twenty-one share a ternary root

Let a finite family of residue classes have distinct odd moduli greater
than one, with \(3,5,7\) absent. Suppose its classes of moduli \(9,15,21\)
agree modulo \(3\). If the family covers the integers, then the union
of the pure \(3\)-power classes has Haar measure greater than \(2/15\).
In particular, \(27\) occurs with a class disjoint from the \(9\)-class.
At least one of \(25,49\) occurs. If \(27\) has the common ternary
root as well, \(25\) occurs.

## Shared ternary root

Use zero clipping at \(5,7\) and the nested upper envelope \(G\) from the
successive restrictions. Retain the \(21\)-charge only on the survivors
of the \(5\)-step. The same argument used to charge the \(35\)-class
outside the \(15\)-class bounds this term by \(w_7\Pr(C_{21}\setminus
C_{15})\) on the product reference law. Write \(w_p=1-b_p\), with
\(b_5=3486/100000\) and \(b_7=81/10000\). The suffix upper envelope,
the two credits and the continuation apply with the caps listed below.
[Successive restriction](84-mixed-class-phase-overlaps.md)

Let \(\lambda\) be the mass of the common root, and let \(s,t\) be the
masses of the quinary root of \(15\) and the septenary root of \(21\).
The probability of \(C_{15}\) is \(\lambda s\), and the probability of
\(C_{21}\setminus C_{15}\) is \(\lambda t(1-s)\). The resulting combined
charge is
\[
\lambda\bigl(w_5 s+w_7t(1-s)\bigr).
\]
It is increasing in \(\lambda,t,s\) on the respective cap intervals,
since \(w_5>w_7t\). If \(s_*,t_*\) are their ordinary depth-one caps,
\(\rho_3\) is the upper density of the ternary reference law, and \(h\)
is an upper bound on the surviving Haar mass of the specified root,
then the joint charge is at most
\[
\rho_3 h\bigl(w_5s_*+w_7t_*(1-s_*)\bigr).
\]
In the complete early comparison, the two separate slots had total
charge \((\rho_3/3)(w_5s_*+w_7t_*)\). Replace that term by the shared
bound, keeping all other slots and the suffix comparison.

## Three cap choices

The \(9\)-class removes a depth-two cylinder from the shared root,
so \(h\le2/9\). If the \(27\)-class is disjoint from the \(9\)-class
and in the same root, use \(h\le 2/9-1/27=5/27\).
A total pure ternary deletion of at most \(2/15\) gives
\(\rho_3=15/13\). Otherwise use \(\rho_3=6/5\).

When \(25\) or \(49\) is absent, the pure deletion at \(p=5\) or
\(p=7\), respectively, is at most \(p^{-2}/(p-1)\).
This gives \(s_*=20/99\) or \(t_*=42/293\), respectively. Without
the corresponding absence use \(s_*=4/19\) and \(t_*=6/41\).
The lower replay for \(b_5\) evaluates the auxiliary envelope after
the \(7\)-step at zero clipping separately for the two septenary caps.
The floor for \(b_7\) uses no earlier coordinates.

| Additional condition | Complete upper numerator at scale \(10^{14}\) | Upper bound with the shared charge |
| --- | ---: | ---: |
| Pure ternary union at most \(2/15\) | \(105051335258083\) | \(3032079236475819623/3038100000000000000\) |
| Both \(25,49\) absent | \(104912804828573\) | \(8671759884987251033/8702100000000000000\) |
| \(25\) absent and a disjoint \(27\)-class shares the root | \(105892617745419\) | \(3623661074857901489/3653100000000000000\) |

Each bound leaves reserve greater than \(1/1000\), which the
[exact verifier](../verify_minimum_nine_modulus_fifteen.py) checks.
The unconditional quadratic continuation gives an uncovered integer
under each additional condition.

With \(27\) absent, the pure ternary union has measure at most
\(3^{-2}+\sum_{e\ge4}3^{-e}=7/54<2/15\). A \(27\)-class contained
in the \(9\)-class may be deleted. This proves the required disjoint
class for a cover and justifies its use in the last implication.
