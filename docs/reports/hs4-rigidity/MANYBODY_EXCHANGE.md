# HS4 many-body exchange: reproducible evidence

The mathematical source is [HIGUCHI_SUDBERY_MANY_BODY_EXCHANGE.md](https://github.com/the-omega-institute/trureturing/blob/lane/theory/hs4-rigidity-chirality-1014/docs/develop/theory/HIGUCHI_SUDBERY_MANY_BODY_EXCHANGE.md) in main PR #16382. It gives written proofs for a two-cell cubic determinant law, a finite Schur-complement error bound, a microscopic entangled ground state and graph-triangle corrections. This program supports those derivations without replacing universal proof by numerical sampling.

## Reproduction

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python verify_manybody_exchange.py > actual.json
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -O verify_manybody_exchange.py > optimized.json
cmp actual.json optimized.json
python -m py_compile verify_manybody_exchange.py
```

Normal and optimized runs produced identical JSON in the reference environment: Python 3.13.5, NumPy 2.3.5, SymPy 1.14.0, SciPy 1.17.0. The actual result is retained in `manybody_exchange_result.json`. The executed source has SHA-256 `20b7b870e9f39d4d38510d757d53611eae8f7565ba321c52b0de213c1b182bbb` and Git blob `09149ff1382cd070d9e4f9f3bf0cff0144722cc4`; the fetched remote source blob agrees. Its numerical output is environment-dependent in the final floating-point digits.

## Exact and numerical scope

The 1,972 exact assertions include every one-spin singlet/triplet relation, all 144 ordered second moments and all 1,728 ordered three-spin Y coefficients. They also include the polynomial cofactor identity in all 16 independent coupling variables, nine complete 256-dimensional physical coupling constructions, the matched-link coefficients, paired-permutation invariance, order-three character projectors and the complete symbolic triangle coefficient. An assertion comparing matrices counts once, not once per matrix entry.

The 65 floating-point assertions use absolute tolerance 2e-10. Seven matched-coupling values g/Delta from 0.001 through 0.08 test the four lowest eigenvalues and energy-dependent Schur complements; eight seeded non-symmetric coupling matrices test the spectral error bound. The three-cell tests use 4,096-dimensional sparse microscopic operators on a path and two weighted cycles. They compare projected second/cubic coefficients, not the full three-cell eigenspectrum. The smallest analytic remainder can be below the floating-point tolerance; numerical passing does not certify that inequality at arbitrary precision. Its universal justification is the resolvent proof.

For g/Delta=0.01 the independently computed full eight-spin spectral gap is approximately 5.075286370733219e-5 Delta and the squared overlap with the encoded Bell vector is 0.9999492479604078. The analytic bounds are respectively 4913/100000000 Delta and 9409/9418. The latter are symbolic conclusions of the proof, not intervals inferred from floating-point eigenvectors.

Three negative controls reject a zero cubic YY coefficient for matched links, a nonzero determinant for the chosen rank-two contrast, and identity action of the three-cycle. They are concrete sanity checks and are not an independent proof-kernel audit or general mutation-testing guarantee. Prior HS4 programs and outputs are unchanged and are not counted as rerun by this program.

## Boundary

Only real isotropic two-spin inter-cell exchange is included in the microscopic Hamiltonian. Second-order states can already have chirality correlations, although the explicit YY coefficient is absent at that order. The cubic determinant governs a specified effective coordinate convention. A singular contrast excludes that cubic coefficient only. Graph triangles have an additional term; edge-wise summation is valid to this order on a forest.

No new Lean compilation, independent mathematical review, full inherited HS4 certificate replay, quantum hardware circuit, DMRG calculation, uniform-volume remainder, phase transition, geometric face-matching equation or factoring speedup is asserted. The Schur transform, Schrieffer-Wolff method and order-three Fourier readout serve different mathematical tasks.
