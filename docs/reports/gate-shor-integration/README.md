# Legal gates, state machines and phase-correct order readout

This is independent evidence for [theory PR #16894](https://github.com/the-omega-institute/trureturing/pull/16894), specifically `docs/develop/theory/FIB_LEGAL_GATES_STATE_MACHINES_AND_SHOR_READOUT.md`. It has no dependency on the HS4 experiments or their PR #11. Earlier HS4 readings are not reused as checks of this line.

## Reproduction and exact source

Run from this directory:

```sh
python verify_gate_shor.py > actual.json
python -O verify_gate_shor.py > optimized.json
cmp actual.json optimized.json
python -m py_compile verify_gate_shor.py
```

Both executions were performed and returned identical JSON in Python 3.13.5, NumPy 2.3.5 and SymPy 1.14.0. `py_compile` succeeded. The actual output is retained in `result.json`.

Executed source: 6,638 bytes, SHA-256 `e73b94e8a6e60cfaa12f57b69efa8da5e040d31aef465e8c506badfad16557e8`, Git blob `1dfd981a1abe6780568e1f1add4387c9ac869bee`. The fetched remote source blob matches the executed local file.

## Coverage

There are 5,033 exact assertions, seven floating-point comparisons and three concrete negative controls. Counts refer to assertions, not matrix entries or a general proof-kernel audit.

The exact corpus includes no11 words of lengths 1 through 8, totaling 141 words and 22 sweep cycles. It checks legal toggle involutions, sweep bijectivity, integer phase exponents modulo four, cycle return-phase consistency and paired phase cancellation. A symbolic three-cycle characteristic polynomial verifies the arbitrary product of its three edge amplitudes. The paired cancellation on this corpus is an integer-phase check; the general identity is proved in the theory.

For the actual five-mode, three-position sweep, the direct matrices satisfy U=-P_tau. Starting from 001 on the classical three-cycle, a six-dimensional Fourier clock has peaks 1,3,5 for U and 0,2,4 for P_tau. Direct simulation of U tensor conjugate(U) on the basis-copy code recovers the latter peaks. The clock computations are floating matrix calculations with absolute tolerance 2e-11. They do not compile or execute a primitive quantum circuit.

The work-register example compares the reduced clock for a clean period-three function with the clock obtained when the entire exponent is retained in an orthogonal work register. The latter becomes maximally mixed and has a uniform Fourier distribution. The universal work-state condition is derived by partial trace in the theory.

Fibonacci modular tests cover every modulus from 2 through 128. Integer arithmetic compares repeated transitions and repeated-squaring powers, checks the faithful (1,0) probe and the nonzero-orbit bound, and verifies primitive periods using the prime factors of these small candidate periods. The test does not assert that factoring arbitrary candidate periods is free or classically efficient. For ideal eigenphase labels, the two-sample denominator-LCM probability is checked exactly for periods 1 through 40 against the Chinese-remainder product.

The three negative controls reject: literal quantum return at the classical three-step period; extension of the paired monomial identity to a Hadamard gate; and preservation of the clean clock by retaining the full exponent as garbage. They use concrete well-separated counterexamples, not a general mutation-testing guarantee.

## Boundaries

The mathematical constructions assume explicitly specified phases, registers and access to required gates or their controlled powers. Tensor-conjugate cancellation is restricted to monomial gates with a supplied conjugate implementation. It is not a construction for an arbitrary unknown unitary or an unknown-state cloning procedure. Generic fusion-chain gates need not have finite order or cheap large powers.

The modular-arithmetic section gives an architecture based on known matrix powers, reversibility and uncomputation. No primitive gate list, optimized resource count, error-corrected realization, QPU experiment, Lean compilation, independent mathematical review, classical-versus-quantum speedup or historical-priority claim is supplied. Universal statements rely on the ordinary proofs and their stated assumptions, not on the finite test count.
