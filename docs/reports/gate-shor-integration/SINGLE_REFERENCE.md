# One phase address and native AURIC certification

Companion to `docs/develop/theory/FIB_LEGAL_GATES_SINGLE_REFERENCE_AND_CERTIFIED_READOUT.md` in main PR #16894. The current dev measurement source is `AURIC_FIB_ATOM_FOUR_CORNER_COHERENCE_AND_SINGLE_COPY_CERTIFICATION.md`, C-117/C-121, read at commit 2c22f027c10d658fbf32942a70f717537c61b2d0. This delivery stays on the independent gate/state-machine line and does not change HS4 experiments.

## Actual reproduction

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python verify_single_reference.py > actual.json
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -O verify_single_reference.py > optimized.json
cmp actual.json optimized.json
python -m py_compile verify_single_reference.py
```

These commands were executed successfully. Normal and optimized output were identical. Reference environment: Python 3.13.5, SymPy 1.14.0, NumPy 2.3.5. The output is retained in `single_reference_result.json`. Floating digits can depend on environment; comparisons use the tolerance below.

Executed source length: 13,014 bytes.
SHA-256: `68159012253ed651001e21c2b1825dd84711b2755ee80f85cb0a9272850170f5`.
Git blob: `dfab8fe9def2d867693563320f0a313519c09628`.

The program emits `single_reference_pulses.json` when run. Its reference SHA-256 is `ff7bbecb30cf95d84f9d16b42b241a454fc286c655a715d18d23df4050a7fdc1`. This generated pulse list is not checked in; it is reproducible from the source.

## Scope of the counts

There are 23 exact symbolic equalities, 48 structural assertions, 55 floating comparisons, and three concrete negative controls. Matrix equalities count once, not once per scalar entry.

Symbolic coverage includes native controlled-Y and SWAP on two and three encoded data qubits, both controlled-Y directions, arbitrary-angle remote rotations for each of two reference locations, all six complete five-mode witness effects, and the exact linear-score/error-shift examples. The full five-mode effects include the apex and bottom-apex coherences; no assumption of an initially bottom-supported state is needed for these identities.

Structural assertions check the seven-pulse CY, twenty-one-pulse SWAP, distance-dependent remote-rotation costs, the witness costs 7 and 28, and enforcement of the unique second-quadrature address. They are finite menu checks, not additional independent proofs of the general compiler theorem.

Floating checks cover every reference/target pair for one through four encoded data qubits, 24 perturbed five-mode measurement comparisons, and a complete five-qubit Fibonacci-modulo-two order probe. For that probe, all 32 input columns are compared against an independently constructed controlled-permutation/inverse-Fourier matrix in the complete 89-dimensional legal space. It uses 26,174 pulses, of which 327 use the second quadrature and all do so at native site five. The operator residual is approximately 5.9681e-14 and final separator leakage norm 2.2557e-14. Matrix-comparison tolerance is 2e-10; the sampled unitary-error/total-variation comparisons allow 1e-12 slack. These are floating checks, not interval certificates or hardware readings.

The three negative controls reject identification of controlled-Y with CNOT, nonpositive null drift after a deliberately biased allowed distributional perturbation without score correction, and the claim that the compiled circuit uses only phase-zero pulses. These are concrete sanity controls, not a proof-kernel audit.

## Mathematical versus operational conclusions

The universal routing and measurement statements rely on explicit operator identities in the theory. The time-uniform probability guarantee is proved with the conditional-null, independent-setting and calibrated error assumptions; it is not established by enumerating finite data streams. A fine-grained final word measurement realizes the stated coarse effects but is not asserted to have the same postmeasurement map as a coarse Lueder instrument.

Only one data address has a second phase quadrature. All native locations still require independently addressed phase-zero control. The minimal-address result concerns the unchanged complex encoding and specified initialization; alternative real-amplitude encodings and supplied reference states have a different contract. Reducing phase-capable addresses increases this unoptimized example's pulse count relative to the previous richer menu.

No new Lean compilation, canonical admission, external independent review, QPU execution, finite-blockade calibration, optimal synthesis, factoring speedup or universal entanglement detection is asserted. Earlier experiment sources and outputs are preserved and are not counted as rerun by this program.
