# Exact checks for the type-D quiver moment-map gap

Run `python3 check.py` with Python 3 and SymPy 1.14.0. The program uses integer,
rational and exact Gaussian-integer algebra. It contains no floating-point
calculation. `check.out` is its standard output; successful execution ends with
`ALL_CHECKS_OK` and exits 0.

The source is John Maar, *On the gap of quiver representations*,
[arXiv:2606.02805v1](https://arxiv.org/abs/2606.02805v1), Proposition 2.1 and
the final paragraph of Section 1.3. Proposition 2.1 is the fidelity check:
the normalized traceless matrix moments satisfy the defining trace identity
for all orientations in ranks 4–6 and three deterministic dimension/matrix
patterns, including zero-dimensional vertices. The source moment map excludes
the zero representation; the program explicitly excludes 16 such cases.

The other checks examine the proposed ordinary proof's ingredients:

- The affine-D Cartan identity and signed-root consequence for every dimension
  vector with coordinates 0, 1 or 2 in ranks 4–8.
- Spanning-tree/unbalanced-edge bases for every family of at most three signed
  roots in ranks 2–4, with additional full root families, zero columns and repeats.
- Exact pseudoinverses and kernel-preserving projections for all 210 unordered
  pairs of nonzero binary affine-D4 dimension vectors of Tits form at most 1.
- Squared moment norms for 510 actual unstable single-arrow representations,
  in affine ranks 4–9 with endpoint dimensions at most 6, compared to the
  proposed universal lower bounds. Rank-deficient endpoint row/column scalings
  put these representations in the SL nullcone.

These finite checks do not prove a universal gap bound, do not search all complex
representations, and do not certify a Lean settlement. The positive settlement
route is a general argument about critical points, scalar moment blocks and signed
incidence matrices. Its intended bound is
`gamma >= 1 / (1765 q^7 d^(11/2))`, hence
`gamma >= 1 / (2000 q^7 d^6)` for affine D, with `d >= 1`.
Finite D follows by giving the affine diagram's leaf 0 dimension zero.

The program and this report are original experiment material. No source TeX or
third-party implementation is reproduced here; SymPy is an external dependency.
