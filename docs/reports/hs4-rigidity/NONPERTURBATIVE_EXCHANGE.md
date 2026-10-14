# HS4 nonperturbative exchange: exact identities and bounded checks

The mathematical source is Sections 11-16 of `docs/develop/theory/HIGUCHI_SUDBERY_MANY_BODY_EXCHANGE.md` in [main PR #16382](https://github.com/the-omega-institute/trureturing/pull/16382). Its new statements concern the full two-cell Hamiltonian, not a finite-order effective Hamiltonian: all-positive matched-link ground-state uniqueness, a uniform-link exact ground block, an ordered Schmidt spectrum, recovery instruments and finite-graph energy brackets.

## Reproduction

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python verify_nonperturbative_exchange.py > actual.json
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -O verify_nonperturbative_exchange.py > optimized.json
cmp actual.json optimized.json
python -m py_compile verify_nonperturbative_exchange.py
```

Normal and optimized executions produced identical JSON in Python 3.13.5, SymPy 1.14.0, NumPy 2.3.5 and SciPy 1.17.0. The actual output is retained in `nonperturbative_exchange_result.json`. The executed source SHA-256 is `ca4f7f8f1549db869eebb66d1cb745941a4c8e5703880497ccb916c345e6815d`, and its Git blob is `fe019ee2f9ff4ab6513f1faf402fbdd56600c066`. The fetched remote source has the same blob. Floating-point final digits can depend on the numerical environment.

## Executed coverage

The 105 exact assertions cover spin reversal, all single-cell spin projectors, the three sandwich-operator identities, the explicit three-dimensional matrix and its characteristic polynomial, the full physical 256-by-3 invariant-subspace action, the weak/strong endpoint identities, and the complete recovery operation on all nine sector matrix units. Recovery checks include Alice's 47-outcome normalization and trace-preserving completions of Bob's conditional maps. These count equality assertions, some of which compare full matrices, rather than individual scalar entries.

The 713 numerical assertions include 37 logarithmically spaced positive ratios g/Delta from 0.001 to 1000, 12 independently disordered positive matched-link cases, and four complete three-cell ground-energy brackets. Each uniform case compares the full 256-dimensional ground eigenvector and energy with the exact three-dimensional formula, computes the entire Schmidt spectrum, and checks the stated ordering and recovery probability. Hamiltonians are scaled to max(Delta,g)=1 before the main floating comparisons. The maximum scaled energy discrepancy is approximately 3.11e-15, and the maximum absolute vector-infidelity residual is approximately 2.67e-15.

For the three-cell path and triangle, the full twelve-spin Hamiltonian is restricted only to total magnetization zero. Its 924-dimensional sector contains a representative of every total-spin multiplet: twelve spin-1/2 sites carry integer total spin and the Hamiltonian is SU(2)-invariant. This permits a complete ground-energy comparison, rather than a singlet-code variational calculation. The sparse eigensolver's returned eigenpair residual is checked. These remain floating computations, not interval-certified eigenvalue bounds.

The ordinary numerical tolerance is 2e-10 and the graph eigenpair residual tolerance is 2e-8. Assertions labelled positive are bounded sanity checks at this tolerance; they do not establish strict positivity over a continuum. The all-parameter uniqueness and majorization statements depend on the written proofs. Three negative controls reject uniqueness at zero coupling, identification of the small block's excitation gap with the full gap, and unit success probability of the strong-coupling spin-zero filter.

## Mathematical and operational boundary

The all-positive unequal-link theorem supplies uniqueness and full Schmidt rank. The exact three-dimensional ground block, majorization in a single coupling parameter and explicit ground-family recovery require uniform matched links. The block contains the exact full ground state but not, in general, the full system's first excited state. Pointwise finite-system uniqueness gives no spectral gap uniform in size or toward zero coupling.

Spin-zero filtering is a flagged, probabilistic protocol. The separate fixed deterministic recovery uses general operations within each four-spin cell and one-way classical communication; some branches change total spin and are outside the exchange-only menu. The matrix checks use fixed output-coordinate representatives of the conditional maps. No physical implementation, gate complexity or hardware experiment is asserted. The deterministic map is guaranteed only on the specified three-dimensional entangled subspace, not on arbitrary eight-spin inputs.

The finite-graph bounds use the exact pair energy through a standard cluster lower bound and a matching-product variational upper bound. They are not an exact many-body phase diagram, a uniform Schrieffer-Wolff error estimate or a thermodynamic gap theorem. No new Lean compilation, full inherited HS4 certificate replay, independent mathematical review, canonical repository admission or historical-priority claim is made. All preceding experimental programs and outputs remain separate prior evidence and are not counted as rerun by this program.
