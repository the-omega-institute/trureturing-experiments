# FIB parent-neighborhood exact data

Volume X §13 proves the existence of an open two-real-parameter neighborhood of `(3/7,-35/29)` and a fixed finite rational b-interval policy. The policy takes the actual sixth-clock M reading, then selects one actual seventh-clock M or J reading and recovers the full original world. No numerical neighborhood radius or interval table is claimed.

The mathematical argument combines the nonnegative compact source envelope, the VI27.3 bound on actual MM witnesses, exclusion of bounded-b pole limits using X2.5, extension of X11/12 to all relevant nonnegative sources, regular inverse coordinates at smooth image points (including F44), and squarefree-fiber b coordinates at the exceptional layers. The Python checks below verify exact polynomial data used by this argument; they do not themselves prove its analytic or policy steps.

All file loading and result output are relative to `__file__`. The directory can be copied as a unit and executed from any working directory. Python 3.12 or newer is required; only its standard library is used, with no dependency installation.

```sh
python3 verify_source_regularity.py
python3 verify_discriminant.py
```

`source_regularity_macaulay.certificate.json` is the exact supplied 47071-byte source certificate, SHA-256 `f4a6f91da3e833ede44549a532a9b32d7695e3e24d744b8e6bac8e5e46f654ab`. The first script reconstructs all five original generators, verifies the multiplier sum equals 1 and maximum product degree 8, and checks the tail/corner coefficients and their nonvanishing on the closed parameter box of radius 1/100. It writes `source_checks.json`.

`coefficient_manifest.json` identifies byte-exact inputs from experiment commit `7c45a21dbea5cdc0974f878a96dc64281cc5db07`, directory `fib_shared_response_176/global_preflight`. The second script checks those local input hashes, the exact retained discriminant factorization through degree 240, factor squarefreeness and pairwise coprimality, and H's constant e-leading coefficient. It reuses the supplied resultant rather than recomputing it and writes `discriminant_checks.json`. The original source back-substitution units, reduced-fiber squarefreeness and separation certificates remain the immutable suppliers cited by X11/12; these scripts do not regenerate them.

The payload has no Lean declarations and supplies no Lean verification, independent review, CI, merge, finite-precision, noise, whole-parameter-plane, alternate-arm, ancestor, original root8/root9, or optimality result. The theorem is ordinary pure mathematics, with the exact-arithmetic checks distinguished from the proof. General local-ring and inverse-function tools are reused literature; the FIB neighborhood continuation is repository-derived.
