# Frame compactifications and MT-algebra realizations

Source: Sebastian D. Melzer, Cerene Rathilal and Ranjitha Raviprakash,
*A choice-free approach to Wallman compactifications*,
[arXiv:2606.05832v1](https://arxiv.org/abs/2606.05832v1), Remark 3.10
(printed PDF page 10; HTML Remark 3.14).

Run with Python 3; the program has no third-party dependencies:

```sh
python3 docs/reports/melzer-rathilal-raviprakash-mt-compactification/check.py
```

`check.out` is the program's output. Subsets are integer bit masks; unions,
intersections, complements and comparisons are exact. No floating point is
used.

The fidelity check reproduces Example 3.8: the Boolean homomorphism from
the two-element algebra to the four-element indiscrete algebra induces an
isomorphism on opens but is not onto on Boolean carriers, so it is not an
MT-compactification.

The realization check includes all 35 topologies on labelled sets of zero
through three points, all 1,225 ordered pairs of these frames, and all
25,577 frame homomorphisms between them. The frame-homomorphism enumerator
is independently checked against all function tables on the 36 pairs
having at most two points in each representative. Each frame element is
checked to be the union of the join-irreducibles below it, which makes the
larger enumeration exhaustive.

For a homomorphism g:K→L, represent K and L as their actual open subsets,
put the graph of g in the Boolean product P(Y)×P(X), and define interior
as the union of all graph elements below its argument. The program checks
the four interior axioms, exact fixed-point set, the source-frame order
isomorphism, the projection's complete Boolean homomorphism laws, the
MT-morphism inequality and the commuting equation. All 904 dense onto
maps yield MT-compactifications, with surjectivity on both carriers and
opens and density on opens. Every open cover is included in the finite
compactness check, including the empty cover of the trivial algebra.

These results support the proposed general graph argument. They do not
prove the unbounded assertion, test nonspatial frames, or verify a general
Funayama representation supplier. The intended positive settlement is a
general proof, not finite enumeration.
