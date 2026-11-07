# Weyl weak identities over finite fields

Source: V. Futorny, P. Koshlukov and J. F. Schwarz, *Identities of Representation
of Lie subalgebras of the first Weyl algebra*, arXiv:2611.06305v1, Section 6,
Lemma 6, equations (8)–(10), Theorem 7 and the finite-field question.

Run with Python 3 and no additional packages:

```sh
python3 check.py > check.out
```

All arithmetic is exact. Weyl multiplication in the integral PBW basis
`P^i Q^j` uses `QP = PQ - 1`; an independent adjacent-word reducer checks the
multiplication rule on all 255 binary words of length at most seven.

The fidelity check reproduces the source's cubic weak identities and the
quartic consequence on all tuples from the basis `{P,Q}` over the integers.
Their multilinearity extends those evaluations to the full generating span.
The same identities are checked exhaustively on the full spans over F₂ and F₃.

The expected refutation is checked for q = 2, 3, 4, 5, 8, 9. For every pair
in the q²-element span, `[u,v]` equals its scalar determinant and
`[u,v]^q - [u,v]` vanishes. In the free associative polynomial, the degree-two
component is the nonzero polynomial `-[x,y]`. Extension fields are constructed
using the irreducible moduli `t²+t+1`, `t³+t+1` over F₂ and `t²+1` over F₃;
the field laws and nonzero inverses are checked exactly.

The scalar-polynomial envelope is checked on every coefficient tuple through
degree five over F₂ and F₃, and degree four over F₄ and F₅: vanishing on all
scalars is equivalent to zero remainder modulo `t^q-t`. Every scalar is
realized as a commutator value.

These finite computations support the uniform argument; they do not prove
membership or nonmembership in the infinite consequence ideal. The latter
uses its general minimum-degree-three filtration. No complete finite-field
basis is claimed. `check.out` is the complete successful output and ends in
`ALL_CHECKS_OK`.
