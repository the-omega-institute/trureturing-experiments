# Exact coefficients for the uncalibrated complex and real PB interface

The mathematical statements, derivations and source conditions belong to
`AURIC_FULL_PB_RAW_CALIBRATION_AND_COMMON_SOURCE_TRANSPORT_CONTINUATION.md`,
sections 56–60, in the main repository. These assets supply finite exact
coefficient data and symbolic identities, rather than physical PB samples.
No result here is Lean kernel certified.

Use Python 3.12 or later and the pinned packages in `requirements.txt`:

```sh
python3 -m venv /path/to/venv
/path/to/venv/bin/python -m pip install -r requirements.txt
/path/to/venv/bin/python complex/normal_form_certificate.py
/path/to/venv/bin/python complex/rational_corner_checks.py
/path/to/venv/bin/python complex/exact_matrix_checks.py
/path/to/venv/bin/python real/real_pb_precheck.py verify
```

All paths resolve relative to the program files. Outputs replace the adjacent
mathematical JSON data. `normal_form_arrays.npz` contains only the two reduced
matrices; their identities are the little-endian uint16 hashes in the JSON,
independent of ZIP container metadata. `real/real_pb_precheck.py build`
regenerates the real finite certificate from the fixed seed; `verify` instead
reconstructs every selected frame and constraint over the rationals and uses
both FLINT and separate plain-Python elimination for the determinant.

The complex map has all 1920 real coefficients, including scalar-input and
scalar-output components. Universal common-corner block preservation has the
152-parameter normal form

$$
D(a,b)=\eta(a,b)I+\ell(a)b-\ell(b)a+i(aHb-bHa)/2.
$$

The 7128 by 1920 matrix from 33 rational Gaussian corners has a 1768-square
minor equal to 551 modulo 1009. Its normal-form model has 152 independent
columns. Rational reconstruction checks every input/readout factor and the
retained arrays. This is a coefficient-space certificate, not a theorem that
33 PB frames or finite physical preparations suffice. Full PB then eliminates
the remaining parameters by the analytic argument in the main text.

The real map has all 450 coefficients. Ninety rational frame/actor pairs give
450 necessary full-kernel rows with determinant 123800 modulo 1000003. The
finite legal zero-map witness has `phi(L)=-10` and budgets `258101/40000`,
`242101/40000`. Mixed spectator kernel vectors and full matrix products are
retained. The proof concerns exactly Sym4(R), O4 and real q, not all real
quantum theories or realifications with an additional phase reference.

The symbolic file also checks two scalar exclusion families, two alternative
families, the corrected compressed identity `4(d1-2d2+d3)/15`, and both free
noncommutative sufficiency identities. The constant product is Z*Z for the
positive orientation and ZZ* for the negative orientation; unstarred ZZ is
not constant. Exact rank, Gram, exterior algebra and modular lifting are
standard tools. The two repositories' existing licenses apply to these
project mathematical programs and data; no external article is redistributed.
