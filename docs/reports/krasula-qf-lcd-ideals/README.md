# LCD ideals in quasi-Frobenius rings

`check.py` is an exact, deterministic finite check for the open question in
Dominik Krasula, *Free and projective LCD codes*, arXiv:2609.15423v1, Section 5.1.
Run it with Python 3:

```sh
python3 docs/reports/krasula-qf-lcd-ideals/check.py
```

The program checks twelve explicitly constructed algebras over the prime
fields F_2 and F_3: truncated polynomial rings, two-by-two matrix rings,
radical-square-zero two-cycle algebras, and three direct products. It checks
multiplication associativity, both distributive laws and the identity against
the exact operation tables; addition is the coordinate vector-space operation. A full-rank matrix for
the associative pairing lambda(xy) certifies each algebra is Frobenius. Every
vector subspace is enumerated once through its reduced row echelon basis;
closure under basis multiplication selects every left and right ideal.

For source fidelity it reproduces Theorem 19, the central-idempotent
classification of LCD ideals of a finite Frobenius ring, and the ideal case of
Theorem 15's size condition. It also checks the proposed annihilator identities,
the contracting double-left-annihilator orbit, and both left and right LCD
classifications on every eligible ideal. The non-symmetric two-cycle example
checks that arbitrary QF ideals need not have equal left and right annihilators
or be fixed by the double left annihilator. A product-factor ideal supplies a
control showing an LCD ideal may contain a nonzero square-zero ideal; a
semiprime-ideal theorem alone would not settle this question.

`check.out` is the program's output. Success is exit 0 with final line
`ALL_CHECKS_OK`. All arithmetic is integer arithmetic modulo a prime; no
floating point, randomness, or third-party Python package is used. These
checks exhaust the ideals of the listed algebras, not all QF rings. They do
not prove the classification for arbitrary or infinite QF rings and do not
supply a Lean proof.

The check program is original project code under this repository's license.
The source is cited for its mathematical statements; no paper text is copied
into the program.
