# HS4 rigidity and parent-Hamiltonian checks

Companion to [trureturing PR #16382](https://github.com/the-omega-institute/trureturing/pull/16382) and `docs/develop/theory/HIGUCHI_SUDBERY_RIGIDITY.md`.

## Reproduction

Run `python verify_hs4_rigidity.py`. The reference environment is Python 3.13.5, SymPy 1.14.0 and NumPy 2.3.5. The script has no network or repository-checkout dependency. It prints JSON and raises an exception on a failed check. Its explicit checks also run under `python -O`; the standard and optimized runs produced identical JSON in the reference environment. Numerical residuals can vary slightly across platforms.

## Evidence

`result.json` records 50 exact checks and 2,240 numerical checks on 560 states: 512 seeded random normalized complex vectors and 48 small perturbations of the Higuchi-Sudbery state. The exact checks include orthogonality of the two conjugates, all six two-body marginal matrices, their convex mixtures, determinant-invariant values, local Pauli actions in Bell-product coordinates, total-spin kernel dimension, scalar-chirality identities, and the complete spectrum of the explicit three-local parent Hamiltonian.

Three false exact identities are deliberately rejected: the reversed chirality eigenvalue, equality of the opposite nonzero LU invariants, and a Gram matrix with conjugation removed. The reported exact-check count includes the check that all three controls are rejected.

The parent spectrum is 0 (multiplicity 1), 2 (4), 3 (3), 4 (3), and 7 (5), in the coefficient normalization specified in the theory. The two labelled-party invariant values are plus/minus sqrt(3)/864.

## Scope

These are independent finite symbolic checks and bounded numerical tests. They do not constitute a Lean replay, a formal proof of the universal Bell-product normal form, or a replacement for the full written classification and stability proofs. The numerical tests use tolerances and do not establish universal inequalities. They do not certify historical priority or physical implementation.

Input constructions are written directly in the script. No external numerical tables or copied proof certificate are required. The original HS4 global certificate is not replayed here.
