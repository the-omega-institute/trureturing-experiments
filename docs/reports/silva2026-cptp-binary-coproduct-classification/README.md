# Full-matrix CPTP binary coproduct classification: exact finite checks

Source: J. M. S. T. da Silva, *On quantum channels: extreme points, topology, and categorical properties*, arXiv:2608.20689v1, Chapter 7 and §7.2. Proposed answer: a binary coproduct exists exactly when at least one input Hilbert space has dimension zero; the other input is then a coproduct.

Run with Python 3.12 or newer and SymPy 1.14:

```sh
python3 check.py
```

The stored output is `output.txt`. All entries are integers or exact rational complex numbers; no floating-point calculation supplies evidence.

The fidelity checks recover the source's hom-set dimension formula on 16 positive dimension pairs by the exact rank of the trace-preservation constraints, with a strictly positive depolarizing Choi matrix witnessing an interior point. They also reproduce the arithmetic obstruction in Corollary 7.3 and enumerate the small Pythagorean pairs left by Theorem 7.2.

The settlement checks use 46 explicit CPTP injection pairs with exact Kraus completeness certificates. They include dimensions (1,1), (3,4) and (5,12), mixed-unitary and pure-preparation maps, rational complex rotations of support subspaces, and additional unused output dimensions. Support projections are calculated from the normalized output matrices rather than assumed. Every input matrix unit is checked for support, discriminator factorization and invariance under relative phase; a cross-support matrix distinguishes relative-phase conjugation from the identity. The zero-input checks cover dimensions 0 through 4.

These finite families validate the identities used in the proposed general argument. They do not enumerate all CPTP maps and do not prove the universal classification. No Lean proof is included.
