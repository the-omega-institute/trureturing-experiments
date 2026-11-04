# C113 finite-error rank certificate

This unit retains exact finite data for the original C113 four-label instrument
and its selected 5 by 5 response table. The complete ordinary mathematical proof
is in C113 sections 37–48 in the main repository. The checker uses exact rational
arithmetic modulo t^8, t=theta^2, and exact SymPy polynomial identities. It uses
no numerical sign samples and supplies no Lean/kernel verification.

## Requirements and source

Python >=3.12 and SymPy 1.14.0 are required. The checker rejects a different
SymPy version. Its required `--repository` is a local clone of
https://github.com/the-omega-institute/trureturing containing commit
`2ecdfbb41b5bd234229ea9b61176a174a4c4de90`. It reads only `git show` at that fixed
commit and path:

`docs/develop/theory/AURIC_INTERNAL_BALANCED_CALIBRATION_AND_FILTER_COMPOSITION_GEOMETRY_C113.md`

The pinned source has 163924 bytes and SHA256
`39921aa9b28ea442cf0e0a20f93a6ff872207baec6738d1cb1cc5f3318e9db5a`.
Both size and hash are checked before calculation. That entire source is also
an unchanged prefix of the C113 source at
`dc97e7044750e3b788268e692f467f00cb2d53c1`, whose 177357 bytes have SHA256
`4da5fcd3dd408293f640b699cfa181f5742933e9aa6005fdf20197cd49b94853`.
The later source includes section 36; its three-design Delta_S contract is
different from the section 32.7 nominal-calibration Delta used in this unit.

## Replay from any working directory

Use absolute paths for the checker and repository. The output destination is
explicit and may be outside both repositories:

```sh
python3 /path/to/trureturing-experiments/experiments/auric/c113-finite-error-rank/certify_rank.py --repository /path/to/trureturing --output /tmp/c113-certificate.json
cmp /tmp/c113-certificate.json /path/to/trureturing-experiments/experiments/auric/c113-finite-error-rank/certificate.json
```

The parent directory of the output file must exist. The checker writes only
that file and prints a compact result. Repository and output paths are not
manifest inputs. The manifest records fixed dependency requirements rather than
host paths, interpreter build strings or timestamps; replay under the stated
requirements reproduces its complete bytes.

## Exact data and mathematical boundary

The manifest retains every characteristic coefficient through t^7, all 25
four-by-four minor coefficient arrays, independent Newton trace polynomials,
all five generalized Vandermonde identities, the Laurent moment/determinant
identity, determinant/compound/singular-value constants, endpoint remainder
ratios, the three small-eigenvalue leading coefficients, matched CP moments,
all nine uniform CP error factors and the Hankel error multiplicities.
The checker asserts the lower-order vanishings, independent Newton consistency,
rational inequalities and all nine CP error factors. Analytic sine, Cauchy and
Taylor remainder arguments and spectral deductions belong to the full owning
source proof, rather than to these finite assertions.

For 0<theta<=2^-13 the original table has
`theta^6/8230118400 <= sigma_min <= (65536/3645) theta^6`, with leading term
`(64/12605) theta^6+O(theta^8)`. Arbitrary complex entry errors strictly below
`theta^6/41150592000` preserve rank five. The rank-four CP comparison is a
hypothetical fixed complete instrument in the assumed complex/Born/CP class;
it supplies no tilted-axis operation to the original menu. Its selected-table
error is below `16 theta^6`, with no order-six full-word error assertion.
The original section 32.7 Delta expression cannot activate the entry threshold;
this compares an upper-bound expression, not the true diamond defect.

No certified statistical accuracy, source stationarity, physical minimal
Hilbert dimension, exact full-future dimension, optimal constants or native
derivation of complex scalars, Born evaluation or CP is supplied. The finite
unit does not settle FourPath or the continuing AURIC goal.
