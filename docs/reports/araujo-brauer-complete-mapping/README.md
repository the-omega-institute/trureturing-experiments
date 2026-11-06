# Brauer and partial Brauer complete mappings

Exact Python checks for the proposed answer to the whole-monoid existence clause of Problem 15.11 in João Araújo, Wolfram Bentz, Peter J. Cameron, Kevin Hendrey and Michael Kinyon, *Complete Mappings of Semigroups*, arXiv:2608.25092v1. For every positive integer n, the proposed answer is that each monoid has a complete mapping exactly when n = 1 or n ≥ 4. Classification of every individual mapping is outside this target.

Run with Python 3, using only the standard library:

```sh
python3 check.py
```

`check.out` is the resulting output; the final line is `ALL_CHECKS_OK`. The program uses integers, finite sets, permutations and exact `fractions.Fraction` arithmetic.

The fidelity check reproduces the source proposition with TeX label `prop:brandt-latin`: complete mappings of the combinatorial Brandt matrix-unit semigroup correspond to Latin squares. An independent exact-cover search finds 2 mappings at order 2 and 12 at order 3 and checks the resulting Latin squares. This Brandt semigroup is distinct from the Brauer monoid.

The settlement checks cover:

- literal full diagram products and explicit complete mappings for B₁, PB₁, B₄ (105 diagrams), and PB₄ (764 diagrams);
- every product in B₂, PB₂, B₃ and PB₃, confirming that a product is a unit exactly when both factors are units, together with exact nonexistence checks for complete mappings of S₂ and S₃;
- telephone-number and perfect-pairing counts through complement size 10, and the rank-2/rank-3 index parity formulas through degree 40;
- the two explicit base pairings of through sets and the four-new-vertices induction through degrees 38 and 39;
- all half-diagram support entries through degree 5, including diagonal identities, inverse transport, and positive balanced weights; degree 6 uses samples of at most 80 actual half diagrams per rank;
- complete anchored permutations at B₆ rank 2 (45 indices, 21 two-cycles and one three-cycle) and B₇ rank 3 (105 indices, 51 two-cycles and one three-cycle), with every cycle equation checked modulo sign.

For the degree-4 constructions, a group complete mapping supplies the good ranks. For bad ranks with even index m, write the half-diagram indices as i, λ ∈ {0,…,m−1}, let ε be the input group element's parity, and put h = floor(2i/m), μ = λ + 2i + ε modulo m. A local bijection from parity-ε elements to parity-h elements is chosen so that the products are distinct. The multiplier has row λ, group coordinate from that local bijection, and column μ. The diagonal sandwich entry is the identity. The program verifies the resulting bijection and product bijection on the actual diagram carrier.

Adjacent through sets guarantee compatibility for Brauer half diagrams. The corresponding assertion for partial Brauer half diagrams fails: in degree 2 the half diagrams with sole through vertices 1 and 2 and singleton complements are incompatible. The partial Brauer proof uses even indices instead and does not need that stronger assertion.

These finite checks support a general combinatorial proof route; they do not prove the unbounded classification, formalize the published lifting theorems, or claim a fresh Lean verification. The complete argument must establish the principal-factor coordinates, all-parameter pairing induction, source lifting criteria, and unit obstruction in Lean after preregistration. No third-party source text or data is redistributed here; the program is original code for this report.
