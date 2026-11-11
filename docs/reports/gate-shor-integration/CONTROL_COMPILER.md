# Native no11 controllability and encoded circuit compilation

This is supporting evidence for Sections 10-15 of `docs/develop/theory/FIB_LEGAL_GATES_STATE_MACHINES_AND_SHOR_READOUT.md` in [main PR #16894](https://github.com/the-omega-institute/trureturing/pull/16894). It is independent of the HS4/many-body experiments. All previous gate-Shor sources and outputs remain unchanged and have not been rerun by this verifier.

## Reproduction and exact executed source

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python verify_control_compiler.py > actual.json
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -O verify_control_compiler.py > optimized.json
cmp actual.json optimized.json
python -m py_compile verify_control_compiler.py
```

These normal and optimized executions returned identical JSON. Python source compilation succeeded. The reference environment is Python 3.13.5, NumPy 2.3.5 and SymPy 1.14.0. Actual output is retained as `control_result.json`. The executed file has 14,144 bytes, SHA-256 `45e0d57b5a43c8d091ad6da2bd178c4b9a7092e5536091260732b2041086cd32`, and Git blob `578f63e2f18c46e5a4818f886a6cbb8e50dcff18`. The fetched remote source blob agrees with the executed local file. Floating residual digits can depend on the execution environment.

## Counts and what they mean

The output's 860 `exact_assertions` comprise **89 symbolic matrix/scalar equality assertions** and **771 exact integer adjacency checks** performed while routing CNOTs. They are not 860 separate symbolic proofs. Symbolic checks cover the encoded isometry, arbitrary-angle single-qubit actions, the separator conditional phase, five-pulse CZ and both nine-pulse CNOT orientations for two and three encoded qubits. They also check the gauged-real generators, antisymmetry, parity conjugation and the orthogonal blind pair for legal lengths 1-5.

There are 458 additional `modular_rank_assertions`, including replay of each stored nested-commutator ancestry, trace-zero checks and rank comparisons at the prime 1009. Integer matrices X_i and K_i=Y_i/i generate the complexified Lie algebra. A nonzero full minor modulo 1009 gives a genuine characteristic-zero rank lower bound for those actual integer Lie words. Trace zero gives the corresponding upper bound. The one-extra-quadrature ranks are 3, 8, 24, 63 and 168 at legal lengths 1-5, exactly d^2-1. All three positions of the phase reference are also tested at length three. The phase-zero ranks observed at those lengths are 1, 3, 10, 28 and 78. The universal phase-zero orthogonal containment and all-n full-control theorem rely on the written proofs and the published open-path toggle-group theorem, not extrapolation from this list.

The six floating comparisons have absolute tolerance 2e-10. Three concrete negative controls check: an orthogonal pair invisible to the phase-zero protocol; failure to reset the separator when replacing its full pi pulse by a pi/2 pulse; and different states with identical weights in one fixed spectral measurement. These controls are not a general proof-kernel audit.

## Actually compiled finite circuit

The source circuit has five ordinary qubits: a three-bit clock and two bits implementing the Fibonacci permutation (x,y)->(x+y,x) modulo two. Controlled powers 1,2,4 are intentionally implemented by repeated small powers, decomposed into Toffoli and then H/phase/CNOT gates. The inverse eight-point Fourier transform is similarly decomposed.

The native compiler produces pulses represented by an explicit position, generator phase and angle. It uses phases 0 and pi/2 at individually addressable sites. Five logical data qubits are embedded on nine physical occupation sites, separated by four initial zero bits. The full legal Hilbert space has dimension 89 and the code dimension is 32. Data-to-data routing is explicit and unoptimized.

The generated source contains 234 H/phase/CNOT gates after the first decomposition. Routing uses 771 nearest-neighbor CNOTs, and the final native sequence contains 7,328 pulses. The verifier compares the complete 89-by-32 encoded action against an independently constructed controlled-permutation/Fourier matrix, on all 32 input columns. Its operator-norm residual is approximately 2.09397e-14 and its final code-leakage norm is approximately 1.36964e-14. It also compares the clock output distribution for initial target 10 and clock 000.

These dense all-input matrix computations are a small-system classical validation and scale exponentially with the number of data qubits. The compiler's pulse-list construction instead has the stated circuit-dependent overhead. Neither claim implies an efficient classical simulator of arbitrary quantum circuits. No circuit optimization, large-integer factorization or hardware run is represented by these data.

## Important scope distinctions

The all-n controllability theorem only needs one site with an extra Y quadrature, alongside all addressed X_i. The explicit short gate compiler uses the supplier's richer menu with Y quadratures at every logical data site. No short method for propagating a single physical phase reference to all compiler gates is claimed.

All operations use the ideal no11-projected pulse model, arbitrary exact angles and no uncontrolled drift. Separators can be excited during a pulse, but every intermediate state remains in the legal Hilbert space and separators return exactly to zero at each compiled gate boundary in the analytical identities. Floating simulation verifies those identities on finite examples without creating an interval or hardware certificate.

Full controllability, informational completeness and short computational synthesis are different conclusions. The general group-control input, the toggle-group theorem, standard CNOT universality and blockade-mediated phase-gate precedent retain their published ownership. The new continuation gives their concrete correspondence and compiler within this declared menu. No new Lean compilation, Scribe, canonical admission, independent mathematical review, historical-priority claim, physical implementation or quantum-complexity advantage is asserted.
