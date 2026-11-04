# Odd commuting-graph invariants

`check.py` checks finite cases of the proposed affirmative answer to Tânia
Paulista, *Commuting graphs of inverse semigroups and completely regular
semigroups*, arXiv:2511.10612v1, Problem 7.1:
<https://arxiv.org/html/2511.10612v1#S7>.

Run with Python 3, using only its standard library:

```sh
python3 check.py
```

All arithmetic uses integers. Prime fields are tabulated exactly; the field
of order 9 is constructed as F₃[t]/(t²−2), with its field laws checked by
exhaustion. Projective matrices are canonical representatives modulo nonzero
scalars. Commutation is tested by direct matrix multiplication modulo scalars.
No numerical approximation is used.

The fidelity check reproduces the constructions in source Theorems 4.3 and
5.3: the proper commuting graph of S₃ × Cₘ has clique and chromatic numbers
2m, for 1 ≤ m ≤ 6. It computes and removes the entire center.

The proposed-family checks cover PGL₂(Fq), q = 3, 5, 7, 9, 11. They verify
all centralizer sizes, the split and nonsplit involution neighbor counts,
the shifted quadratic character sum, the three-stage deletion ordering, an
explicit q-coloring, and a nonsplit-torus clique of size q. A clique and a
proper coloring with the same cardinality certify both invariants exactly.
For q = 3, 5, 7 and m = 3, 5, the program independently recomputes every
commutation pair in PGL₂(Fq) × Cₘ, removes the entire product center, and
checks the blow-up adjacency formula and matching clique/coloring witnesses.

`check.out` is the program's complete output. Its last line is
`ALL_CHECKS_OK`. These are finite checks of a proposed general argument;
they do not prove the universal assertion, and no Lean proof is included.
The program is original code; the cited source supplies the problem and
fidelity target, not program text.
