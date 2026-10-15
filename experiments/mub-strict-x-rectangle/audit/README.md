# Independent exact audit

Run from the experiment directory:

```sh
python3 audit/verify_independent.py \
  --result results/delta-2neg24 \
  --input input.json \
  --preconditioners audit/preconditioners.txt \
  --output audit/independent_verification.json
```

Only Python 3's standard library is required. Run without `-O`. The source files
are expected beside `input.json`; `--source DIRECTORY` can override that location.
The result directory must include its 32 chart reports, combined report, matrix
bounds, root enclosures, and `centre_preflight/centers.txt`.

The checker:

- verifies source/data hashes and the minimal root-radius change;
- audits all 32 completed chart reports and closed-tree counts;
- recomputes all 60 enlarged root contraction and strict-inclusion bounds using
  independent arbitrary-precision integer interval arithmetic;
- performs three independent root-enclosure refinements with fixed dyadic
  preconditioners and checks all 1,770 distinct-ray comparisons and 1,656
  forbidden pairs on those independently derived enclosures;
- additionally checks the C++-reported sharpened enclosures;
- verifies the exact center Hadamard Grams and shift commutation over
  Q(i,sqrt(21));
- freshly validates the center roots and canonical shift fixation, using the
  explicit center matrix rather than an old success report.

`preconditioners.txt` contains proposed dyadic matrices with denominator 2^40.
They are untrusted inputs: no inverse accuracy is assumed, and every use is
validated by fresh exact contraction and inclusion bounds. Altering a proposal
cannot create mathematical acceptance unless the exact inequalities still hold.

The connected-root-branch propagation lemma is required to transfer the freshly
verified center fixation across the enlarged rectangle. See the experiment's
proof documentation. This audit does not independently replay a full global
cover trace; it checks fresh cover reports against the reviewed source and
independently validates the local roots, center symmetry, and graph margins.
It is not Lean-kernel verification and makes no novelty claim.

The JSON contains hashes rather than machine-specific filesystem paths, so its
provenance is stable when the experiment directory is relocated.
