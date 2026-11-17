# Fixed-parent shared-response certificates

This bundle reproduces the ordinary mathematical certificates in
`FIB_SELF_CALIBRATING_RELATION_MATHEMATICS_X.md` sections 9--10. The parent is
`(theta, xi) = (3/7, -35/29)`, with literal `MMJMM` readings
`(0, 0, 0, 6, 72/7, 2700/203)`. The model retains every strictly positive,
rank-one, trace-one source, every nonzero gain polynomial of degree at most
three, arbitrary constant drift, the original clock, immediate readings and
one trajectory that cannot be reset or copied. The mathematical dependencies
are VI sections 8, 23, 24 and X section 8 at
[`064f028cd737ae0241e2042aa34e5d99c8669707`](https://github.com/the-omega-institute/trureturing/tree/064f028cd737ae0241e2042aa34e5d99c8669707/docs/develop/theory).
The X sections 9--10 source input has SHA-256
`b8b44bd1922cf610b5c4c18a6c7443fd523caff7c8c81be6705b24d70b245d6d`.

Python 3.12 or newer, SymPy 1.14.0 and python-flint 0.9.0 are required;
`requirements.txt` pins the Python packages. No network access is used by the
reproduction. From the experiment repository root, with these packages available:

```sh
python3 fib_shared_response_176/reproduce.py /absolute/path/to/FIB_SELF_CALIBRATING_RELATION_MATHEMATICS_X.md
```

The command runs the six entrypoints below in a temporary directory inside the
bundle, checks the advertised data hashes, compares regenerated polynomial bytes
and rational isolations, and leaves the preserved coefficient data untouched.
It prints six program rows with `"exit": 0`, followed by:

```json
{"passed": true, "global_identity_checks": 17, "bprime_checks": 78, "necessary_candidates": 112}
```

`--output /absolute/path/to/result.json` also writes compact mathematical result
data; `results.json` contains the verified values and Bprime interval measurements.
The standalone suite programs write coefficient and result files beside their
own source; use `reproduce.py` to keep the advertised data bytes unchanged.

| Entrypoint | Exact output and scope |
| --- | --- |
| `global_preflight/probe.py` | Constructs `R`, `E`, the chart polynomials, exceptional resultants and the image resultant. `R` has degrees `(u,b)=(16,4)` and 85 terms; `E` has degrees `(u,b,e)=(12,3,1)`; the single `e`-dependent factor `H_factor_1` has degrees `(b,e)=(16,16)` and 153 terms. |
| `global_preflight/discriminant.py` | Constructs `Res_e(H,H_e)`, a nonzero polynomial of degree 240 in `b`, and its factorization. |
| `global_preflight/support.py` | Constructs the monic squarefree collision support of degree 202, including leading-coefficient, content, chart and discriminant exceptions. Factor degrees are `1,2,2,2,4,24,32,44,91`; exact Sturm counts are `1,2,2,2,4,20,24,8,49`. |
| `global_preflight/isolate_support.py` | Isolates all 112 distinct real support roots by rational intervals of width at most `10^-12`. Regenerated intervals and multiplicities agree with `support_isolations.json`; elapsed-time fields are excluded from that comparison. |
| `global_preflight/verify_certificate.py` | Checks 17 exact identities and support properties: the child line, `Res_n(Q_A,T_A)=A^4 R`, the chart denominator, cancellation, image resultant, factorization, discriminant, squarefreeness and exceptional-factor coverage. |
| `bprime_display_correction/verify_bprime.py` | Checks 78 identities and literal-display/interval conditions over the entire Bprime band `[16423371109/10^10,1642337111/10^9]`. There are exactly two legal complete worlds in each such M-child fibre and their actual terminal MJ response intervals are disjoint. Arithmetic uses symbolic rationals and outward integer intervals at scale `10^80`. |

All five advertised files are preserved byte-for-byte in `global_preflight`:

| File | SHA-256 |
| --- | --- |
| `R.json` | `a7b5c02552dea6e43bbd4919bb0f0db88e92e838fca26b1bbc976603437d568e` |
| `E.json` | `de87cab361a10a0f7a84a15a3738b5e8cad7702331698b6fcf2f3b1e530aff36` |
| `H_factor_1.json` | `a4c611369cd5eedc73c7314581a507da880b7db30f5d2fc6fd01eb9866d5229e` |
| `collision_support.json` | `ce0e91ac8fcf7b8637f0520f368938a9f165d305cc76c156c949cce5bbb727f1` |
| `support_isolations.json` | `9f0cd27bde73dcbedac70eca973ee9942ee7d8bc46998232ec168e471ab62cd7` |

The other retained coefficient files are inputs to the certificate verifier,
support construction and isolation: `A`, `D`, `cancelled`, `phi_reconstructed`,
`H_raw`, `H_main_discriminant`, the three `res_R_*` exceptions and nine
`support_factor_*` polynomials. Coefficient JSON uses variables `(u,b,e)`;
each term contains an exponent tuple and an exact rational coefficient string.
The large raw discriminant and intermediate factor dumps are reproducible
temporary outputs and are not inputs to the retained entrypoints.

`source_provenance.json` records the project-owned source paths and original
SHA-256 values. These input sources are not contained in experiment base
`3963df3a954ef46316d4ddcecd735ada7b754ef0`; no source commit is claimed for them.
`scan_fibers.py` supplies the symbolic matrices, resultant and chart expressions
imported by `probe.py`; its optional numerical scan is not run.
`dependencies/fib_x_section8_display_correction/verify_section8.py` supplies the
exact interval arithmetic and mathematical-display parsers imported by the Bprime
verifier. Its source bytes are unchanged. The two importing programs resolve
these dependencies relative to the bundle; all other copied programs are
byte-identical to their provenance inputs. `reproduce.py` only orchestrates these
programs, compares results and preserves the advertised bytes.

The 112 roots are only necessary candidate responses for MM collisions. No
physical specialization at all 112 candidates, emptiness of
`C_MM intersection C_MJ`, root eight/nine-read decision or nine-read upper bound
is established here. At an exceptional candidate, the original positive-real
equations and all exceptional branches still require specialization. Each
terminal action may use its own pair of worlds, but the two pairs must share the
same actual sixth response within this M arm. No cross-arm equality of sixth
responses is imposed. This executable reproduction is by the implementation
worker in the same model family; it is neither independent review nor Lean
verification.
