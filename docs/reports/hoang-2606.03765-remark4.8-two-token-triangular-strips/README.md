# Two-token triangular strips

Exact checks for Duc A. Hoang, *Token-sliding realizability for complements,
Cartesian-products, and grid graph families*, arXiv:2606.03765v1, Remark 4.8
and Table 2. The proposed settlement is that no finite simple graph has
its full independent-two-set token-sliding graph isomorphic to
C₃ □ Pₙ when n ≥ 6.

Run with Python 3.12 or newer:

```sh
python3 check.py
```

The standard-library program uses integers, sets and exact backtracking.
Its stored output is `check.out`; successful completion ends in
`ALL_CHECKS_OK`.

- Source fidelity: the hosts K₃ ⊔ K₁ and K₃ ⊔ K₂ realize the first two
  strips, as in Corollary 4.4. The searches for n = 3, 4, 5 reproduce the
  exclusions in Lemma 4.6.
- Literal independent-set/sliding adjacency agrees with the complement-root
  edge implementation on every labelled host of orders 0–5.
- Brute enumeration covers all 33,868 labelled hosts of orders 0–6. It finds
  84 realizations of the first strip and 70 of the second, counted over these
  labelled hosts and orders; none realizes strips 3–5.
- A separate search assigns distinct root edges to target vertices, keeps
  required and forbidden root edges, and uses a necessary bipartite matching
  test for pending root edges. It uses neither the proposed induced-cycle
  obstruction nor the claw/interface classification. Its decisions on all
  ten connected target isomorphism classes through four vertices agree with
  exhaustive host enumeration through order five.
- For connected targets with m states, active root support is connected:
  target adjacencies share root endpoints. Its m edges therefore use at most
  m+1 vertices. Isolated root vertices correspond to universal host vertices
  and can be discarded. Fresh endpoint names are introduced canonically, so
  the root-label search covers every host order up to this irrelevant padding.
  For n = 6 and 7 this gives support bounds 19 and 22, respectively. The
  searches visit 65 and 97 nodes and find no realization.
- Direct chordless-cycle enumeration finds only the n triangles and 3(n−1)
  squares in the first ten strips. A separate check verifies the induced
  six-cycle on the six blockers forced by a same-center interface; a C₆
  positive control checks that the cycle enumerator detects six-cycles.

These finite computations are evidence for the proposed structural argument,
not a Lean proof of the uniform statement. No Lean proof is included here.
The source's existing local interface argument and a new induced-cycle
obstruction are the proposed route to that uniform statement.
