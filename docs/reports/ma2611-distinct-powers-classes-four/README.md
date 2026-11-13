# Distinct powers with generator classes of size four

Exact finite checks for Dingying Ma, *On Distinct Powers in Finite Groups*,
arXiv:2611.13582v1, Theorem 1.1, Remark 2.2 and Section 9.
The expected settlement is the affirmative size-four analogue, uniformly in
finite groups and positive n. These computations do not prove that assertion.

Run with Python 3.12 or newer, using only the standard library:

```sh
python3 check.py
```

`check.out` is the complete output. Success ends in `ALL_CHECKS_OK`.

- Fidelity: all 24 labelled placements in Remark 2.2 fail; Lemma 4.1 has
  exactly the six stated exceptional multisets over F7.
- Latin arrays: exhaust all restricted-growth-labelled prefixes of orders
  three and four, then all possible final rows blocking every transversal.
  There are no obstructions of order three and 594 labelled obstructions of
  order four, in exactly the two trisotopy classes displayed by Best,
  Hendrey, Wanless, Wilson and Wood, Section 4, arXiv:1612.09443v1;
  published DOI 10.1002/jcd.21566. No external program was reused.
- Four-block assignments: exhaust all four-element subsets of the unit
  group for admissible moduli m ≤ 31 and m = 49, all positive r ≤ 11 with
  r + 4 below the least prime factor. For m = 121 the unit domain is
  explicitly {1,2,3,4,12,23,34,45}. Reductions modulo a prime may coincide.
- Group checks: exhaust every eligible n-subset, n ≤ 6, in every odd cyclic
  group of order at most 31, C3 × C3, C5 × C5, C9 × C3, C3 × C3 × C3,
  C7 semidirect C3, C9 semidirect C3 and the order-27 Heisenberg group.
  Powers in additive-coordinate models mean repeated group multiplication.
  Every test group's associativity, identities and inverses are checked.
- To exercise larger n with genuine size-four classes, exhaust all nonempty
  subsets of the nine-point domain {0} ∪ {(u,0),(0,u) : 1 ≤ u ≤ 4} in
  C11 × C11. This checks 511 subsets through n = 9, not all subsets of that
  121-element group.

All evidence uses integer arithmetic and `fractions.Fraction`; there is no
floating point or random sampling. Assignment witnesses are checked again
against the power table. The program and this report are independently
authored for the repository; mathematical sources are cited above.
