# Reflection-product counts in floretion triangular coordinates

Exact checks for Creighton Dement, *Equilateral completion in floretion triangular coordinates: locality, product points, and reflection symmetry*, [arXiv:2608.15925v5](https://arxiv.org/abs/2608.15925v5), Remark 7.8 and Table 3. Source version dated 29 October 2026.

Run with Python 3 and the standard library:

```sh
python3 docs/reports/dement-2608.15925v5-remark7.8-all-order-reflection-product-counts/check.py
```

The program evaluates the source centroid sum with integers. It scans every fixed vertex and every exchanged word pair for each of the three reflections, and compares all three squared distances. It tests locality by the synchronized cyclic orbit and product-centeredness by the unsigned quaternion product. All six columns of Table 3 are reproduced for orders 1 through 6. A separate scan of every unordered triple cross-checks the reflection scan, total triangle count and local count at orders 1 through 3.

For orders 2 through 6 the two unknown columns agree with

- T(n) = 8 * 3**(n-1) - 5 * 2**n + 2;
- U(n) = 8 * 3**(n-1) - 3 * F(2*n+2).

At order 1 the counts are T(1) = 1 and U(1) = 0.

The program additionally checks all 64 inverse-transition substitutions for the proposed eight-state certificate, using exact rational arithmetic. Exactly 24 integral parity-compatible edges survive. Exact matrix multiplication verifies M**3 * w = 3 * M**2 * w. Unpruned integer-state dynamic programming independently matches the accepted path counts through order 10.

These checks support a proposed uniform proof; they do not constitute a Lean proof or settle the all-order claim. The closed forms, axis reduction and bounded-carry approach are prior art in Korea Superintelligence Labs, [*Product-centered reflection-symmetric equilateral triangles in floretion triangular coordinates*](https://lean.up.railway.app/papers/floretion-counts-v1.pdf), Conjecture 4.5 and Sections 3 and 6. That paper reports formal finite checks through order 10 and nonformal checks through order 12; its all-order automaton correspondence is explicitly unproved. The proposed additional obligation is a complete backward-closure correspondence followed by induction for arbitrary lengths.

`output.txt` contains the complete successful result, ending in ALL_CHECKS_OK. SHA256 of `check.py`: a2dca55abe6f832eb0b3b6d9b5eedd648c12c4f867b087f449dba6abc42779ad.
