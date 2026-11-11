# FIB common-base incidence: exact local obstruction

The original problem asks for one entire child with acquired readings
`(0,24,147/5,24,18,p,v)`, `p != 24`, containing independently selected,
within-pair distinct equal-J and equal-M pairs at literal clock 7. Both base
coordinates remain free. The finite discovery sample supplies no admissible
certified four-endpoint solution and no global exclusion.

The substantive exact result is Section 9 of
`FIB_STRICT_REAL_UNEQUAL_DIAGONAL_ACTUAL_M_PAIR.md` in the main repository:
the entire child at its existing exact M-pair base contains exactly those two
worlds. Their actual J answers differ. More generally, throughout the explicit
radius-`10^-24` rectangle of Definition 9.1, every entire child has at most two
worlds and its actual terminal-J response is injective. This is a local
exclusion of the independent-pair predicate over both base coordinates.

The proof retains the existing unsaturated factorization
`Res_d(f,K) = 256 D^8 F0^2 R` and the linear pseudo-remainder in the drift.
Uniform rational contractions cover all 18 complex R roots and all three F0
roots. The positive candidates have nonzero inverse denominators; strict
source tests leave only R slots 12 and 13, with disjoint J ranges. Each slot
determines at most one whole world. This inventory fixes no terminal-M value
and imposes no source-width bound. Gain degree drops, zero samples, repeated
gain roots, zero drift, singular sources and exceptional source layers enter
the necessity argument. Their exclusion where applicable follows from the
certificate, rather than a domain restriction.

| Artifact | Mathematical role |
| --- | --- |
| `paper_certificate.json` | Exact rational data for the printed Section 9 certificate, including all centres, uniform root bounds, inverse filters and J separation. |
| `verify_paper_certificate.py` | Rational-only replay without a root finder; also compares printed centres, source filters and bounds with the source. |
| `existing_M_recheck.json`, `verify_exact.py` | Recheck the existing M pair inside its original defining box and its uniqueness there; independently verify all original literal clock 0–7 equations for one whole cubic-or-lower gain and fixed drift. |
| `existing_M_base_entire_child.json`, `certify_M_base_child.py` | Higher-precision root-centre supplier and entire-child certificate construction. Numerical roots choose centres; rational inequalities certify the result. |
| `inverse_coefficients.json`, `reuse_inverse.py` | Exact reused R, alpha and zeta coefficient normalization and symbolic identities. |
| `negative_source_joint_intersection.json` | Exact three-endpoint joint intersection at a base near `(8.30114690304,76.37500654862)`. Endpoint 0 has `r < 0` and `b < 0`; this is inadmissible. |
| `disjoint_guard_failure_intersection.json` | Exact disjoint joint intersection at a base near `(24.60337353395,133.91984484175)`. M endpoints 2 and 3 have respectively `r > 1` and `r < 0`; this is inadmissible. |
| `certify_intersections.py`, `exact_intervals.py` | Certificate construction and outward integer-grid rational intervals with forward derivatives. |
| `preregistration.json` | Full target, finite workload envelope and falsifiers. |
| `curves.json`, `parent_curves.json`, `pair_catalogue.json` | Numerical pair curves and candidate source pairs for further searches. These do not certify complete components or child inventories. |
| `joint_candidates.json`, `joint_catalogue_candidates.json`, `joint_tb_candidates.json`, `overlap_candidates.json`, `same_pair_candidates.json` | Numerical candidate and guard-failure data from 12,000 total major starts. Small residuals do not certify existence. |
| `companion_data.json`, `positive_child_scan.json` | Numerical companion data for selected bases. These are not exhaustive inventories. |

The reusable discovery programs accompanying these data are `discover.py`,
`full_parent.py`, `trace_parent.py`, `joint_search.py`, `pair_catalogue.py`,
`joint_catalogue_search.py`, `overlap_search.py`, `joint_tb_search.py`,
`same_pair_search.py`, `companions.py` and `positive_child_scan.py`.

The source circuit and symbolic carrier machinery are retained under
`../fib-self-calibration/` at their original owning relative paths. Every
consumer resolves them from its own `__file__`; the current working directory
and the names of sibling checkouts do not select suppliers. No mathematical
supplier has been moved, renamed or reimplemented. The three retained source
files are byte-for-byte original copies. The [dependency inventory](DEPENDENCIES.json)
records every retained file's original and publication SHA-256, original
checkout and revision, consumer edges, and licensing. All 35 original
experiment files are retained; all original JSON files remain byte-identical.
Modified program headers identify the publication path/input changes.

The exact replays need Python >= 3.9, SymPy and mpmath. Tested versions are
CPython 3.9.6, SymPy 1.14.0 and mpmath 1.3.0, pinned in `requirements.txt`.
The retained numerical discovery programs additionally need NumPy 2.0.2 and
SciPy 1.13.1, pinned in `requirements-discovery.txt`; the exact replay does
not import NumPy or SciPy. These are package dependencies, not vendored source.
The original project files use the repository's Apache-2.0 [LICENSE](../../LICENSE);
the supplier checkout has the identical license and no NOTICE file.
SymPy, mpmath, NumPy and SciPy distribute their own BSD-style license notices
with their packages.

Run from the private experiment repository root, using an interpreter with
`requirements.txt` installed (for example,
`python3 -m pip install -r experiments/fib-common-base-6ada6767-attempt-1/requirements.txt`).
Provide the public theorem as an explicit read-only input:

```sh
THEORY=/path/to/main/docs/develop/theory/FIB_STRICT_REAL_UNEQUAL_DIAGONAL_ACTUAL_M_PAIR.md
OPENBLAS_NUM_THREADS=1 python3 -B experiments/fib-common-base-6ada6767-attempt-1/replay.py --theory "$THEORY"
```

The approved theorem input SHA-256 is
`f39a62ba731197f79ba97288739282ba69af7f4602e4b18692beb9f8c309253d`.
Both source-comparison entry points require this hash and fail on a different
input. The theorem itself is not copied into this private experiment unit.
The `replay.py` wrapper verifies the inventory before and after replay,
disables bytecode output, and installs a Python audit hook that rejects file
access under both original experiment checkouts listed in the inventory. It
runs the three exact entry points below, replays the retained negative-J seed,
and resolves the discovery-only quantitative-arc `blocks` input without
starting any discovery program. Its terminal portability summary lists the
three supplier files actually read inside this checkout. The wrapper adds no
mathematical implementation.

The equivalent standalone exact entry points are:

```sh
OPENBLAS_NUM_THREADS=1 python3 -B experiments/fib-common-base-6ada6767-attempt-1/reuse_inverse.py
OPENBLAS_NUM_THREADS=1 python3 -B experiments/fib-common-base-6ada6767-attempt-1/verify_exact.py --theory "$THEORY"
OPENBLAS_NUM_THREADS=1 python3 -B experiments/fib-common-base-6ada6767-attempt-1/verify_paper_certificate.py --theory "$THEORY"
```

Keep Python assertions enabled: do not use `-O` or `PYTHONOPTIMIZE`.
`reuse_inverse.py` regenerates `inverse_coefficients.json` from the original
symbolic supplier. `verify_paper_certificate.py` regenerates
`paper_certificate.json` from the stored rational centres. These outputs must
have the unchanged hashes in the inventory. The other original JSON results
are read-only inputs to replay. The wrapper checks all retained output and
data hashes after execution.

Current exact replay results are the carrier factorization
`Res_d(f,K) = 256 D^8 F0^2 R` (512 terms; degrees `(18,6,6)`), the identity
`A^3 K = alpha*d + zeta + f*Q` (280 alpha terms and 405 zeta terms), all three
stored rational-box contractions, and the original literal matrix/gain
identities at clocks 0 through 7. The sorted-term R normalization SHA-256 is
`ff6b65b00b5d251082949e2318c2a53fb3f583a86c94b3d22c3e15c9ebe038e2`.
The printed certificate replay covers 18 R roots and 3 F0 roots, retains only
R slots 12 and 13 as strict possible sources, checks the printed centres and
source bounds, and gives
`-3076889/100000 <= J_13-J_12 <= -384611/12500 < 0`.
Its leading-R enclosure is
`[1873044541936855184250,1873044541936855184251]`, with the same nonzero
inverse-pole margins recorded in `paper_certificate.json`. The two rejected
joint intersections retain their original source-guard failures.

The replays use exact rational and symbolic arithmetic, without numerical
root solvers. The retained construction programs (`certify_M_base_child.py`
and `certify_intersections.py`) use root solvers to propose centres; they are
not needed to reproduce the exact certificates and are not run for a layout
repair. Discovery scripts retain their original finite workload and seed
settings; several execute searches at top level, so importing them is not a
replay command. The 12,000 major starts remain bounded, nonexhaustive numerical
data. No new numerical search is part of the publication checks. These are
ordinary mathematical and rational-computation checks, with no Lean
certification claim.

The remaining gap is either an admissible exact joint intersection outside
the excluded rectangle, including the possibility of other source components,
or an exclusion covering all original base coordinates and source strata.
The finite search, inadmissible intersections and local inventory establish
neither alternative. The continuing common-base research owner retains that
problem, with both p and v free and overlapping pairs allowed.
