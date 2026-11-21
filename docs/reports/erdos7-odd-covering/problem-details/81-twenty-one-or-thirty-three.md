# A cover with no modulus below nine requires twenty-one or thirty-three

Every finite family of residue classes with distinct odd moduli greater
than one that omits \(3,5,7,21,33\) leaves an integer uncovered. In
particular, a covering of least modulus \(9\) without modulus \(21\)
must have a class of modulus \(33\). If modulus \(11\) occurs, its class
cannot contain the class of modulus \(33\), since such a class of modulus
\(33\) could be deleted without changing the union.

A least-modulus-\(9\) cover omitting \(21\) must also contain \(25\)
and a \(27\)-class disjoint from its \(9\)-class. [Pure-power constraints](82-pure-prime-power-constraints.md)

Use the clipping schedule, the prime-\(5\) charge coefficient \(1-b\)
with \(b=3/100\), and the quadratic continuation from the preceding
proof. With \(21\) and \(33\) absent, the only new omission in the suffix
comparison is the depth-one slot with cofactor \(3\) at prime \(11\).
In the scalar comparison this subtracts \((10/11)/d_{11}\) from the
completed mixed weight when the auxiliary ternary depth is positive.
In the state \(D=1\) immediately after \(5\), the ternary depth is zero.
It remains zero throughout the suffix as the divisor count increases.
Consequently the same lower replay through \(211\) gives the pointwise
floor \(b\) for the upper suffix envelope.
[Prime-\(5\) comparison and floor](80-aligned-nine-and-fifteen-force-21.md)

Temporarily omit the \(15\)-slot at \(5\) as in that proof. The upper
replay applies both omissions \(21,33\) in the positive ternary branch
and gives an unweighted joint upper value
\[
U_0=\frac{91995799802337}{10^{14}} .
\]
The \(5\)-stage charge of the other slots in that replay is still
\(7/95\). The \(15\)-class, if present, contributes at most
\((2/5)(4/19)=8/95\), for any original residue. Thus the cover charge
plus terminal moment on survivors is at most
\[
U_0-\frac3{100}\frac7{95}+\frac{97}{100}\frac8{95}
<\frac{1999}{2000}.
\]
The unconditional quadratic continuation then gives an uncovered
integer. The [exact verifier](../verify_minimum_nine_modulus_fifteen.py)
checks the stated reserve and the earlier comparisons.
