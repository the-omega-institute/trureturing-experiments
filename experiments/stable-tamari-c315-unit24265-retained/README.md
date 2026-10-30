# Stable Tamari Conjecture 3.15: finite experimental evidence

Status: unresolved. The closed probe reports 208,362 distinct directly enumerated
vectors and 1,715,878 product-derived readouts with no nonunimodal entire `F_b`.
These are finite observations, not a proof, a global counterexample exclusion,
an independently replayed result, a kernel certificate, or a delivered resolution.
Retention checked source-file hashes and data structure only; it ran no experiment.

## Fixed source target

Herman Chau, Spencer Daugherty, Sean Grate, Mark D. Haiman, Juan Carlos Martínez
Mori, and Anna Pun, *The Stable Tamari Lattice*,
[arXiv:2610.39257v1](https://arxiv.org/abs/2610.39257v1), Conjecture 3.15.
The paper is CC BY 4.0; see [source.json](source.json) for attribution, clause
locators, recorded source hashes and provenance limits. No paper files are vendored.

For positive `n` and `a,b` in nonnegative `N^n`, define `i <=_a j` when `i=j`,
or `i<j` and `a_k>a_j` for every `i<=k<j`. Define `a <=_T b` by coordinatewise
`a_i<=b_i` and `a_i-a_j<=b_i-b_j` for every `i <=_a j`, using integer differences.
For comparable `p,q`, `ht(p,q)` is the longest chain length, counting edges.
The target is `F_b(t)=sum_{p<=_T q<=_T b} t^ht(p,q)`: the complete coefficient
list, including degree zero, must have a weak peak, allowing ties.
Vectors are neither sorted nor restricted to positive coordinates. The auxiliary
`f_b(t)=sum_{p<=_T b}t^ht(p,b)` is a different polynomial.

## Retained contents and result scope

- `probe.cpp`: literal comparator, source forest comparator/height readout,
  complete finite-box or ideal enumeration, entire `F` accumulation and weak peaks.
  `box` and `checked` inspect every endpoint vector in the box; `single` and
  `singlechecked` restrict to the full ideal of the supplied endpoint. `checked`
  variants compare source forest and literal order; `control` also checks heights
  by longest-chain dynamic programming. `dump` writes endpoint/coefficient rows.
- `independent_reference.py`: separately implemented literal order and longest-chain
  dynamic programming for five small boxes; invokes `./probe` for comparison.
  Both implementations came from the same solver, not an independent evaluator.
- `product.cpp`: exact convolution and weak-unimodality checks using 128-bit counts.
  `product_check.py`: assembles the exact deduplicated product input from the three
  coefficient tables and targeted results. These use only standard libraries.
- `box-results.json` and `extended-box-results.json`: aggregate exhaustive-box
  readouts plus the final endpoint's full coefficients, not every endpoint's list.
  Taking the larger bound per dimension gives `(n,K)` = `(1,80), (2,100), (3,30),
  (4,12), (5,6), (6,4), (7,3), (8,2), (9,2), (10,1), (11,1), (12,1), (13,1),
  (14,1), (15,1)`, with every vector in `[0,K]^n`: 208,206 distinct endpoints.
- `targeted-checked-results.json`: full lists for 240 targeted endpoints,
  including 156 outside the exhaustive boxes. Seed 315261039257 is retained as
  dataset provenance. `last_coefficients` is the entire `F`; `top_f_coefficients`
  is auxiliary `f`, never a replacement for `F`.
- `independent-controls.json`: five literal-order/longest-chain control summaries,
  totaling 233 complete-polynomial comparisons. `dp_height_controls=0` in other
  result files does not assert a dynamic-programming height check there.
- `polynomials3.txt`, `polynomials4.txt`, `polynomials5.txt`: complete coefficient
  tables for `[0,8]^3`, `[0,5]^4`, `[0,3]^5`; each row is comma-separated endpoint,
  semicolon, comma-separated coefficients in ascending degree (trailing commas).
- `product-input.txt`: the 1,852 deduplicated polynomials actually used, with
  representative endpoints; retained to preserve the precise convolution dataset.
  First line is record count; each record is coefficient count, coefficients,
  endpoint dimension, endpoint coordinates. `product-results.json` records no
  failure in 1,715,878 unordered pairs with repetition. Source Corollary 3.14
  identifies a product with `F_(b,0,c)`; these are not direct poset enumerations.

The source forest/height and product identities are documentary premises, not
kernel-checked here. The supplied probe reports that `f_(8,8,4,4)` fails weak
unimodality while its entire `F` passes. Its reproduced source example
`F_(3,2,1)=[14,21,15,10,5,2,1]` also defeats a universal log-concavity route.
Neither observation solves the full conjecture. A uniform whole-`F` weak-peak
argument or an exact full-`F` finite failure, with independent replay and faithful
kernel closure, is still missing.

## Reproduction interface

No reproduction was executed for retention. Work in a separate disposable directory
containing copies of this subtree so programs cannot overwrite retained data.
Use a C++17 compiler and Python 3.13 (standard libraries only):

```sh
c++ -O3 -std=c++17 probe.cpp -o probe
c++ -O3 -std=c++17 product.cpp -o product
./probe singlechecked 3 2 1
python3 independent_reference.py
python3 product_check.py
./product < product-input.txt
```

`./probe checked K K ...` reproduces the recorded box mode; `./probe singlechecked
b1 b2 ...` reproduces a targeted endpoint. No timestamp-bound historical batch
wrapper is needed. Large bounds require a separate resource decision: the program
materializes the coordinate box before its 35,000-vertex guard, and fixed-width
arithmetic is not a general arbitrary-input guarantee. The reported finite runs
used exact 64-bit counts within their stated bounds; convolution uses 128-bit counts.

## Integrity and licensing

`retention-manifest.json` gives source and retained SHA-256 hashes and exact data
projection rules. Four JSON files omit process fields `mode`, `elapsed_seconds`,
and `excluded`; all other values and array ordering are preserved. The original
`excluded` list was empty. Formatting changed; mathematical data were not
independently validated by that projection. Remaining retained source/data are
byte-identical to the named inventory. Original scratch files were not modified.

Original probe programs are Apache-2.0 under the repository's [LICENSE](../../LICENSE),
according to the supplied documentary record. The paper's CC BY 4.0 attribution
is separate. The authors' GPL-3.0 `orders.py` was only compared after independent
implementation; no GPL code was copied, executed, or included, per that record.
