# Shared-response PB rigidity certificate

Exact finite linear algebra and symbolic central-residual witnesses for the full
four-dimensional PB model. Python 3.12 or newer is required. From the repository
root, install `experiments/auric-shared-d-frame-central-pb/requirements.txt` in the
chosen environment, then run:

```sh
python3 experiments/auric-shared-d-frame-central-pb/frame_kernel.py 48
python3 experiments/auric-shared-d-frame-central-pb/central_pb.py
python3 experiments/auric-shared-d-frame-central-pb-precheck/retain_certificate.py
```

The first two sources construct the modular matrix and derive the real-Gram
witness identities. They are retained byte for byte. The separate
`../auric-shared-d-frame-central-pb-precheck/independent_precheck.py` reconstructs
the frames over rational real/imaginary matrices, checks adjoints before
reduction, checks the tensor ordering by direct full-matrix evaluation, and
derives all eight orientation/actor/phase witness cases by real Cauchy–Binet.
It imports neither candidate program. Both implementations use FLINT and SymPy;
this is executable reconstruction with shared backends, not separate proof
kernels or a Lean certificate.

`retain_certificate.py` runs that reconstruction and serializes the mathematical
arrays and `certificate.json` in the precheck directory's `data/`. It overwrites
only these deterministic result files. Their hashes, source hashes, dimensions,
prime, chosen minor and witness polynomials are in `certificate.json`.

All binary arrays contain unsigned little-endian 16-bit integers in row-major
order. The 2304 by 1575 matrix has rank 1575 modulo 1009. Its maximal minor uses
all columns and the 1575 increasing zero-based row indices in `pivot-rows.u16le`;
its determinant is 796. `denominators.u16le` lists the 192 Householder norms.
`input-wedges.u16le` is the 192 by 105 wedge matrix, of rank 105. The rational
construction and ordering are also specified in the main theory, sections 48–50
of `AURIC_FULL_PB_RAW_CALIBRATION_AND_COMMON_SOURCE_TRANSPORT_CONTINUATION.md`.

The matrix eliminates only the traceless-output quotient. The full linear
constraint kernel has dimension 105 and consists of central-valued two-forms;
the symbolic PB witnesses are essential to exclude it. Forty-eight mathematical
frames do not mean forty-eight measurements: every frame retains all nonzero
two-leg complex inputs, all actor probes and all positive full-algebra
functionals. Hashes identify bytes; the matrix construction, exact minor and
rational lift carry the finite argument. No minimum-frame, stability, physical
availability, historical priority or full quantum-origin claim is made.
