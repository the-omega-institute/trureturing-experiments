# Octahedral incidence-ideal boundary dimension

`check.py` supplies exact evidence for the proposed negative answer to the
boundary-component sentence of Betti–Grate–Holleben–Salizzoni Question 3.21,
*Incidence toric ideals and three-point functions*, arXiv:2605.23849v1.

Run with Python 3 and SymPy:

```sh
python3 docs/reports/betti-grate-holleben-salizzoni-octahedral/check.py
```

The output is stored in `check-output.txt`. Every numerical calculation uses
integers, rational fractions, or exact SymPy matrices.

The source checks reproduce its displayed incidence matrix for (4,3,2), full
ranks 15 and 56 for (6,3,2) and (8,4,3), and equality of the exponent vectors
of its displayed quartic and sextic. Generator checks verify parity-half
supports and equal incidence, including coned hyperoctahedra. A separate
64-outcome example checks the alteration calculation and deletion operation.

For k=4, t=3, n=a^7 and p=1/(4a^4), the exact alteration lower bound
p*binom(n,4)-n^8*p^8 exceeds binom(n,3) at a=4,5,8,16. At a=4 it is
2,926,662,418,428 against 732,873,539,584. A coordinate support avoiding both
parity-half monomials of every generator would therefore yield a projective
coordinate subspace of dimension at least 2,926,662,418,427, larger than the
toric dimension 732,873,539,583.

This is evidence for a finite-probability existence argument. The large
support is neither enumerated nor stored. The general-k constant and exponent
checks run for k=4 through 12; these finite checks do not prove the universal
inequalities. No Lean proof, associated-prime classification, projective-degree
calculation, or settlement of the separate degree-equality question is claimed.
