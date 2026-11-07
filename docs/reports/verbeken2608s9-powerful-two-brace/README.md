# Powerful multiplicative 2-groups and nonterminating brace right series

The candidate is the final question in §9 of Brecht F. Verbeken,
“Powerful multiplicative groups do not force right nilpotence in finite braces”,
[arXiv:2608.01884v2](https://arxiv.org/abs/2608.01884v2).

Run `python3 check.py`. The program uses only Python integer arithmetic and the
standard library; its exact output is in `check.out`.

The source fidelity check reproduces its local-algebra and principal-unit
identities at the odd prime 3: dimension 8, ideal order 2187, fourth ideal power
spanned by the nonzero top element, fifth ideal power zero, and nine cubes in
the principal-unit group.

The proposed characteristic-two family has
B = F2[x,y]/(xy, x^(m+1)+y^(m+1)), J = (x,y), m = 4r, and
D(x) = y^m, D(y) = x^m. The invariant character is computed by triangular
principal-unit factorization. For r = 1 the program checks all 512 elements
and all 262144 pairs, including affine composition on an additive basis,
the explicit coefficient formula for the character, inverses, commutators,
fourth powers, and the right series. It obtains G' = G^4 of order four,
class two, exponent eight, and right-series dimensions 9,7,5,3,3.
The nonzero three-dimensional space T satisfies T star J = T, supplying the
finite invariant needed for the mathematical induction at every series index.

For r = 2,3 the program checks the multilinear algebra identities on a basis
and 4096 deterministic integer-generated pairs per parameter. These samples
support the family proposal; they do not prove its unbounded quantifier.
For r > 1 only G' <= G^4 is proposed; G' = G^4 is not claimed.
The general construction, character existence and right-series induction must
be proved by the subsequent formalization. This report contains no Lean proof.

The program is independently authored from the paper's mathematical formulas.
No third-party program or Lean source is copied. Repository licensing applies.
