# Tucker group and the ordinary cube

This independent exact audit concerns Oliver Knill, *Remarks about the Moebius-Kantor graph*, arXiv:2605.30799v1, Sections 1.7 and 3.1. Run with Python 3:

```sh
python3 docs/reports/knill2605-tucker-rubik/check.py
```

`check.py` is original code using only the Python standard library and integer arithmetic. `check.out` is its complete deterministic output. Vertices 0–7 are u₀–u₇ and 8–15 are v₀–v₇. Exhaustive distance-preserving bijection search visits 1,249 nodes and finds all 96 graph automorphisms, reproducing the source's order statement. Every automorphism commutes with the unique distance-four partner map. Its signed action on the eight partner pairs gives 96 distinct cubie images satisfying the twist, flip and permutation-parity constraints. All 9,216 multiplication pairs agree. A generating set of the automorphism group is computed and every generator image passes the constraints.

The six face moves are independently compared to exact integer rotations of physical corner and edge stickers. The general signed-permutation map is additionally checked exhaustively for 0–3 pairs (2,373 multiplication pairs).

These are finite checks, not a Lean settlement. The interpretation of legal coordinates as reachable cube states uses the classical sufficiency theorem; this program does not supply face-turn words for the embedding images or prove that theorem. The intended proof uses the intrinsic partner map and a general embedding for at most eleven pairs, rather than positive enumeration of the 96 automorphisms. The cube has fixed centres and ignores their orientation. Multiplication means function composition, with the right factor acting first.
