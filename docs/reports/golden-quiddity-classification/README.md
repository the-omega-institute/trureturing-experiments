# Golden quiddity classification: exact finite checks

`check.py` checks the unnumbered conjecture after Proposition 3.9 of Flavien
Mabilat, *λ-quiddité et sous-groupes engendrés par un nombre algébrique*,
arXiv:2212.03142v3. The additive coefficient set is φℤ, with φ²=φ+1.

Run with Python 3 and its standard library:

```sh
python3 docs/reports/golden-quiddity-classification/check.py
```

`output.txt` contains the deterministic output. Success requires exit 0 and the
final line `ALL_CHECKS_OK`. No floating point arithmetic or external dependency
is used. Pairs of integers represent the basis 1, φ exactly.

The fidelity checks reproduce the three splice examples after Definition 2.2,
the two signed golden pentagons in Remarks (ii), and Proposition 3.1 on all
integer tuples in [-3,3] of lengths 1–4. Additional identity checks cover U⁵=-I
and zero contraction, and check the zero quadruples for 43 signed parameters,
including 81-digit parameters.

The classification check exhausts all 2,441,405 tuples in {-2,-1,0,1,2} of
lengths 1–9. It finds 1,740 scalar-monodromy words; every one is independently
recomputed as a 4×4 integer matrix using the regular representation. Exactly
nine quadruples and the two signed pentagons are irreducible in this window.
Every quiddity of length at least six in the window admits the proposed zero
contraction or signed-pentagon splice, with both shorter factors checked.

Irreducibility is decided independently of that reduction: for every cyclic
rotation and reversal and every possible internal factor length, the program
solves both endpoints of a hypothetical splice factor. If its interior product
is [[p,q],[r,s]], scalar monodromy forces p=±1, x=-q/p, and y=r/p. These tests
recover every possible factor with entries in φℤ; endpoints are not truncated
to the enumeration alphabet.

These are finite checks of a proposed proof, not a proof of the unbounded
conjecture. No Lean proof or disproof is part of this report.
