# Full-input 48-frame certificate

With Python 3.12 or newer, install `requirements.txt` and run from the repository root:

```sh
python3 experiments/auric-ambient-pb-corner-research/nonunital_frame_certificate.py
```

The adapter uses the relative repository file
`experiments/auric-shared-d-frame-central-pb/frame_kernel.py`, whose immutable
source is commit `651b00b5abcfc10de377ae180bcfde38752bdffd` and SHA256 is
`8f6224dc5c64959c57d375f82f10f6f6605eec11232b86c0614ae9261b6d0331`.
It checks these source bytes before execution; no other checkout or Git history
is required. This Apache-2.0 repository retains the original source and license.
The adapter changes only the input coordinates and exterior-square index range.
Frames, output basis, prime and pivot selection retain their original ordering.

`nonunital-certificate.json` specifies the four little-endian uint16 arrays,
their shapes, hashes, ordering and rational lift. The restriction matrix has
shape 2304 by 1800 and rank 1800 modulo 1009; its specified maximal minor is 20.
The input-wedge matrix has shape 192 by 120 and rank 120. The program regenerates
these canonical arrays and certificate. Exact matrix definitions and their
conditional PB application belong to
`AURIC_FULL_PB_RAW_CALIBRATION_AND_COMMON_SOURCE_TRANSPORT_CONTINUATION.md`
§§53–54 in the main repository.

The separate `../auric-ambient-pb-corner-precheck/independent_precheck.py`
reconstructs the matrix without executing the adapter or its source.
Both programs use numpy/python-flint; independence of construction is not
independence of numerical backend. These finite calculations are not Lean
certification, universal PB sampling, or native-source reconstruction.
