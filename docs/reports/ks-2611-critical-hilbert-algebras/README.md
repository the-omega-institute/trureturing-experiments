# Critical finite Hilbert algebras with zero

The exact check compares the proposed criticality criterion for arXiv:2611.05928v1 with the definition of criticality. It also reproduces the cardinality, height and quasiequation W′₂ of the source's family S⁽ⁿ⁾ (Proposition 4.1).

Run with Python 3 and no dependencies:

```sh
python3 check.py
```

The program enumerates all posets on at most six points, up to isomorphism, by enumerating transitive relations compatible with a linear extension and canonically relabelling them. For each poset it enumerates every hereditary family of antichains containing all singletons and the full set of maximal points. Families related by a poset automorphism may occur more than once. These are counts of representations, not distinct algebra isomorphism classes. The representation theorem implies that all algebras with at most seven elements occur: their representing poset has at most six points, since its singletons and the empty antichain are distinct algebra elements.

The direct test enumerates maps from the implication table, preserving both constants. Every proper image is verified to be a subalgebra, and its map onto that image is a possible coordinate of a product embedding. If the selected coordinates jointly separate every pair, their product is an embedding into proper subalgebras. A verdict of criticality exhausts all endomorphisms. This test does not use the proposed fold criterion. On algebras with at most five elements the homomorphism enumeration is independently compared with all possible functions.

The criterion test separately enumerates idempotent maps on the maximal points and checks the inclusion-maximal faces. Additional exact tests check the idempotent-power normalization of the proper witnesses found in the least-nonmaximal-point case. All arithmetic and comparisons are exact integers, finite sets and equalities. No floating-point computation supplies mathematical evidence.

The completed census contains 70,839 representations over 406 posets with zero through six points; 1,160 representations are critical. All direct verdicts agree with the criterion. The program also checks 6,022 idempotent normalizations of proper witnesses and tests S⁽ⁿ⁾ for 2 ≤ n ≤ 8. Its last line must be `ALL_CHECKS_OK`. The source does not itself assert criticality of S⁽ⁿ⁾; that is a separate checked application of the proposed criterion. The finite checks do not prove the general classification. The full positive settlement requires the general fold-normalization argument.

Source: Tomasz Kowalski and Katarzyna Słomczyńska, *Rule extensions of the implication–negation fragment of intuitionistic logic*, [arXiv:2611.05928v1](https://arxiv.org/abs/2611.05928v1), Sections 2.1–2.3, 3.2 and 4.2.
