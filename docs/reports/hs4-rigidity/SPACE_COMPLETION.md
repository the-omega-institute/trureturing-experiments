# HS4 spatial completion: reproducible finite evidence

The theory source is Sections 9-13 of `docs/develop/theory/HIGUCHI_SUDBERY_AURIC_NETWORK_CLOSURE.md` in [main PR #16382](https://github.com/the-omega-institute/trureturing/pull/16382). The proofs concern an exactly specified spin-1/2 singlet cell, its boundary observation map, a logical-qubit channel family, and finite tensor products of those cells.

## Reproduction

From this directory run:

```sh
python verify_space_completion.py > actual.json
python -O verify_space_completion.py > optimized.json
cmp actual.json optimized.json
python -m py_compile verify_space_completion.py
```

The reference environment is Python 3.13.5, SymPy 1.14.0 and NumPy 2.3.5. Normal and optimized runs returned identical JSON, reproduced in `space_completion_result.json`; compilation of the Python source succeeded. The seed is fixed in the program. Exact-check counts refer to equality assertions, some of which compare complete matrices, rather than to individual matrix entries.

Executed source SHA-256: `f92e96bcb69b81d55535b4541e95977ae02343a7db554384b529d44bca3a22c4`.
Executed source Git blob: `e34ad5d1224c139a5f835d1b40bf8829ba8e5490`.

## Coverage

The program executes 36 exact symbolic checks. These cover the singlet isometry and spin closure, all six compressed pair swaps, all six complete pair marginals for a symbolic logical Bloch vector, signed-volume compression and its squared spectrum, the positivity ball and boundary disk, a shape commutator, binary-entropy curvature, the full symbolic Choi characteristic polynomial, an exact negative Choi witness, an attaining completely positive approximation, the labelled classical POVM, and the two-node cat marginals.

The 5,200 numerical checks use 400 logical-cell states, half pure and half mixed, and 200 three-node states. The latter consist of 100 random rank-at-most-four density matrices and 100 perturbed HS-product states with a small maximally mixed component. All 100 near-product cases give a nonzero product-fidelity lower bound. Matrix comparisons and inequalities use absolute tolerance `2e-11`; entropy sums omit probabilities at most `1e-15`. These are explicit bounded tests, not universal numerical certification.

Three negative controls reject a zero signed-volume operator, the nonnegative-Choi claim for exact orientation erasure, and purity of the maximally mixed logical cell. They are sanity controls for those concrete false assertions, not a general mutation-testing guarantee.

## Scope

The signed triple-product eigenvalues are opposite while the reduced unsigned volume is scalar on this cell. The sharper square-root stability estimate requires exact singlet support and pure input. The retention bound of one half concerns the specified same-qubit unital symmetric channel family with exact erasure; classical labelled observations obey a different contract. The node-product bound keeps local impurity separate from average pair entropy. A spin-network contraction, independent-source preparation, smooth-spacetime limit and gravitational equation are not implemented.

No Lean compilation, inherited HS4 certificate replay, canonical main-repository admission, independent mathematical review, historical-priority claim or hardware realization is asserted. The universal statements depend on the written mathematical arguments, not on the number of finite tests. All earlier experiment sources and outputs are separate preserved evidence and are not counted as rerun by this program.
