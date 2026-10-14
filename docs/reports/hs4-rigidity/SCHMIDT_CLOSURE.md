# Direct Schmidt-block classification and all-state stability

Mathematical source: [HIGUCHI_SUDBERY_RIGIDITY.md, Sections 8-10](https://github.com/the-omega-institute/trureturing/blob/2adffa66e403cf8dbc972745e7831c27b9cd5f78/docs/develop/theory/HIGUCHI_SUDBERY_RIGIDITY.md). Main continuation: [trureturing PR #16382](https://github.com/the-omega-institute/trureturing/pull/16382).

Run `python verify_schmidt_closure.py` from this directory, or pass its complete path from another working directory. Dependencies: NumPy and SymPy. The reference environment is Python 3.13.5, NumPy 2.3.5 and SymPy 1.14.0. The standard and `python -O` executions produced identical JSON; `python -m py_compile verify_schmidt_closure.py` passed. Decimal residuals may vary with a different numerical linear algebra backend.

The reference output is `schmidt_closure_result.json`. This is separate from the preceding `verify_hs4_rigidity.py` and `result.json`; no old checker or old result is overwritten.

## Scope

Twelve exact SymPy checks cover orthonormal and spin-flip identities of the explicitly chosen magic basis; both canonical HS representatives; the collective quaternion rotation action, orthogonality and determinant as polynomial identities; the realness-defect polynomial; the full spin bilinear for a symbolic coefficient matrix; and opposite determinant-invariant values for the conjugate representatives.

The 3,601 numerical checks consist of 640 checks on 160 orthogonal-block states from both O(3) components with product-local rotations; 400 checks on 200 U(3) triplet blocks and their real orthogonal roundings; 2,001 scalar gap samples; and 560 checks on 140 small-deficit states. The small-deficit states are selected by the stated theorem threshold from 210 seeded perturbations. No claim is made about discarded larger-deficit samples. Exact counts are produced by the checker.

Three controls detect the invalid omission of the spin constraint, the invalid omission of the single-cut spectrum condition, and the invalid equality of the two labelled-party determinant invariants. The last control concerns a genuine LU-invariant distinction, not merely orthogonality of two vectors.

The comparator rotations inside the stability construction may act jointly on a qubit pair. They construct a nearby auxiliary vector for an error estimate; they are not asserted to be product-local state-preparation operations. The actual endpoint is in one of the two LU orbits by the separate written classification proof.

These are finite and symbolic computer-algebra checks. They do not compile Lean, certify universal spectral or singular-value arguments, replace an independent mathematical review, replay the inherited HS4 degree-six certificate, or establish historical novelty.
