# Positive inertia of the shifted C27 four-generator core

For the fixed 1263-mode source defined below, one outward dyadic interval
elimination at eta=1 completes all 1263 singleton pivots with strictly
positive lower endpoints. The resulting inertia is

\[
\operatorname{In}(C-\epsilon I)=(1263,0,0),\qquad
\epsilon=27852372618143874844/2^{512}+2^{-542}.
\]

This is a finite numerical certificate conditional on the inherited
coherent-source realization and displacement identities. Authentication of
the saved arrays does not prove those identities or physical realization.
It is not a Lean theorem or official acceptance. It supplies neither a
downstream eta budget nor the feature/paid inverse actions, final bordered
inertia, retained/cofinal positivity, Robin's inequality, or RH.

## Fixed source

The original nonnegative labels are those in the existing good-mode metric
packet. Their order is preserved in `plan.json` and `inertia.json`; no modes
are discarded. The original target is the compressed blob
`2245e47ffd0397a781885a82c6f48585bdd54802` at experiments commit
`f9610d3a5a732abaf68e733e50a7e34e66148ced`, path
`docs/reports/weil-c27/joint-gamma-6765-rho32-v2-target.json.zlib`.
Its SHA256 is
`7546a915368b4a40e0e41ea0ea693403c239598c7a6e96b7b8b07e67e22e98c6`.
The decoded metric SHA256 is
`63bbcebbe00e69b64bf301654a7369ee8c7bb7ecf21dc589f9b408925f8f0d49`.
The common integer-source identifier is
`1ecb02691d9b99e4c181c3f8955ec84dd41328f277bb59bdffcf6d80412dd55e`.

Write I for the saved ordered good set. With the original arithmetic
Schur core H, actual finite middle residual R, and its diagonal weights W,
the coherent source satisfies

\[
C=(H-R^TWR)_{II}-(a_{\rm orig}+\tfrac1{32}a_{\rm reuse})I.
\]

The two padding occurrences are distinct: `a_orig` is the original
subtracted allowance, and `a_reuse` is the upper allowance used by the
joint-boundary consumer. They are not identified or removed. The complete
outer epsilon above is an additional shift. The exact same C is recovered
from the saved diagonal ell and all 131 feature columns by

\[
C_{nn}=\ell_n+(F K_tF^T)_{nn}.
\]

All 112 far and 19 boundary features participate in this recovery. The
off-diagonal source uses its first four saved generator columns G and
original squared nodes:

\[
C_{ij}=\frac{G_i\Omega_4G_j^T}{i^2-j^2},\qquad
\Omega_4=\sigma\begin{pmatrix}J&-J\\-J&0\end{pmatrix},\quad
J=\begin{pmatrix}0&1\\-1&0\end{pmatrix},\quad\sigma=2/\pi.
\]

This is displacement rank at most four, not matrix rank at most four.
The residual Gram enters with a minus sign. These identities are the
explicit source contract; they are not inferred from independently chosen
interval endpoints. The target producers are
`structured_inverse_residual_interval.py` and `joint_gamma_boundary.py`
at the fixed experiments commit. The saved target/metric producer
identifiers and original-input bindings are checked by `interval_core.py`.

## Exact interval computation

`signed_interval_replay.py` authenticates the fixed target and decoded
metric, validates their dimensions, source/basis metadata, full feature
kernel and mode partition, then follows exactly the supplied covering plan.
Eta is the rational 1 and the arithmetic grid is 1024 bits. Input
enclosures remain on their original 640-bit grid before exact lifting;
lifting does not improve the original source information.

For an exact current Schur block A and cross rows B, the independent
diagonal and generator updates are

\[
d'_i=d_i-B_iA^{-1}B_i^T,\qquad
G'_i=G_i-B_iA^{-1}G_K.
\]

The same skew kernel is retained. Signed 1x1 pivots and nonsingular 2x2
pivots are supported. Negative scalar pivots continue normally. Every
reciprocal interval must strictly exclude zero. A 2x2 block uses a strict
determinant sign; positive determinant additionally uses a strict trace
sign. An uncertain block gives no complete inertia. The exact block
identity and inclusion argument are given in `structure.md`.

All elementary endpoints use integer outward rounding. Multiplication uses
corner extrema; squaring handles zero crossing; signed reciprocal uses
the reversed endpoint order. The complete shared 131-feature quadratic
is summed before its single final outward rounding. No floating
eigenvalue estimate, unverified numerical inverse, target producer, or replacement source is
used to establish the recorded signs.

The completed computation took 196.5770605 seconds in its measured
process, with native macOS maximum RSS 599195648 bytes. The actual process
exit and its host carrier exit were both zero. These are observations of
this execution, not estimates or resource guarantees. `inertia.json`
contains the exact intervals for every pivot and the input bindings. Its
SHA256 is
`2a0f2e748634c355d24624e5dbc4b1f87da3afca4d5651894c235cb7665544ac`.

The smallest recorded pivot lower endpoint is approximately
`6.455462e-27`. This is a pivot diagnostic, **not a lower eigenvalue bound**;
congruence does not preserve Euclidean eigenvalues.

## Reproduction

Use Python 3.9 or later. Supply the fixed compressed target and the decoded
metric bytes with the hashes above. If the saved metric is compressed,
decode that single existing packet first; this is not a producer replay.
Run from the repository root, choosing an output path that does not yet
exist:

```sh
python3 docs/reports/weil-c27/signed-core-eta1/run_signed_replay.py \
  --core docs/reports/weil-c27/signed-core-eta1/interval_core.py \
  --target-compressed docs/reports/weil-c27/joint-gamma-6765-rho32-v2-target.json.zlib \
  --metric-decoded /path/to/full-gamma-finite-good-block-metric.json \
  --plan docs/reports/weil-c27/signed-core-eta1/plan.json \
  --output /path/to/new-inertia.json
```

`interval_core.py` is retained byte-for-byte for its authenticated
arithmetic and source decoder. Its separate diagonal-query CLI is not
part of this experiment. The replay CLI requires all input paths and
selects no eta or pivot plan. The saved plan SHA256 is
`dc010d3b6ff069d67bf168786a275806233d49adc7be74710329cf945019b8b3`.
The returned used-plan hash, eta, grid, full ordered coverage, parsed
terminal JSON and actual process exit must all match before accepting a
reproduction. Exclusive output creation prevents overwrites but is not
atomic publication of a completed JSON file.

Source review and small exact dense comparisons exercised signed
reciprocals, negative pivots and zero-diagonal 2x2 blocks before this
application. The complete transcript was separately checked for all 1263
steps, strict pivot signs, source/plan bindings and actual terminal exits.
Those checks do not establish the inherited arithmetic common-source
premises, and this package makes no such claim.

The executed replay source SHA256 is
`8cdfce25b0342e0b0f5b208f186aef1633e69209383a178e42311b76c689256e`.
The packaged replay has an updated module description; all non-docstring
syntax is identical, and its loader pins the packaged source hash. The
saved plan's structure-proof identifier names the inherited conditional
source contract. The program compares that identifier but does not read
or kernel-check a mathematical proof. `structure.md` gives the algebraic
contract and inclusion proof for this package.
