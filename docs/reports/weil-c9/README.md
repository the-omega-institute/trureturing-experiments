# Actual c9 retained comparison and source-moment verification

This package supplies scientific inputs for one sufficient Weil comparison at the actual support half-width $a=\log3$. It reuses Liu's published analytic formulas and source-only moment oracle, described in the [Library note](https://github.com/the-omega-institute/trureturing/blob/dev/Library/Weil/liu2026tailcompensation.md). The author implementation remains an external, attributed dependency. No author retained matrix or old positivity certificate is reproduced here; this package claims neither a new moment algorithm nor a proof of RH.

The project moment packet comes from executing the complete published 6400-bit oracle, then enclosing its 1,024 outputs on the $2^{-1024}$ grid. The affected 256-mode kernel was regenerated at 3072 bits. Every one of its 65,536 exact matrix centers, its tail-column centers, all error allowances and the complete sufficient target agree with the earlier comparison. Its previously recorded positive sign is therefore reused by exact target identity. **No fresh sign execution is asserted for the new manifest.** Reuse fixes the prior certified dyadic trial; it does not assert that a new floating proposal reproduces those trial integers.

[`source-input-result.json`](source-input-result.json) records the source calculation, separate-session verification and target identity. [`result-summary.json`](result-summary.json) remains the original sign evidence and retains its original manifest binding. Independent implementation, independent prime-matrix regeneration, Lean verification and cofinal support positivity remain absent.

## Replay

Use CPython **3.13.12**, `python-flint==0.9.0` and `numpy==2.5.3`. These checked version strings do not attest the interpreter or arithmetic-library binaries. From the repository root, with `python` selecting that interpreter:

```sh
python -m venv build/weil-c9-venv
build/weil-c9-venv/bin/python -m pip install -r docs/reports/weil-c9/requirements.txt
curl --fail --location https://raw.githubusercontent.com/luciferyu666/certified-weil-positivity/b6cd2183c1e79c6c27a34267812a7b2d73ed1b59/frozen-source/reproduction/source/01-check.py --output build/liu-01-check.py
build/weil-c9-venv/bin/python docs/reports/weil-c9/reproduce.py --moment-checker build/liu-01-check.py --out-dir build/weil-c9-reproduction
```

The external checker is pinned to SHA-256 `4ed7661f1f869ed6d37d7bc52f1a0889a8d305a790663dc0ae43d76f2ea0e356`. Obtain it subject to the publisher's [rights statement](https://github.com/luciferyu666/certified-weil-positivity/blob/b6cd2183c1e79c6c27a34267812a7b2d73ed1b59/RIGHTS.md), which authorizes independent reproduction without granting a blanket redistribution license. Author source, manuscript and old certificates are not included here.

The wrapper verifies the scientific program hashes and packed and decoded input hashes. It decodes all three inputs, executes the complete external source-only moment oracle, checks equality of the fresh scientific rows with the retained project rows, and runs the directed consumer at 1536 bits. The oracle's four built-in negative controls must reject their mutations. Use a fresh output directory: the source wrapper refuses to overwrite an existing moment result.

For bindings alone:

```sh
build/weil-c9-venv/bin/python docs/reports/weil-c9/reproduce.py --check-inputs-only --out-dir build/weil-c9-bindings
```

That mode performs neither moment generation nor a sign calculation. A full user-invoked replay writes `result.json`, including the exact dyadic trial and all directed pivot enclosures; it is a fresh sign calculation only if that execution actually succeeds. A successful sign run exits zero. The scientific sign condition is **all 256 directed pivots strictly positive and no directed negative witness**. Elapsed times and floating proposals are not byte-reproducibility requirements.

The original separate native Codex CLI session executed the same finite-sign program on the original hash-bound inputs. It checked all 256 serialized pivot lower endpoints and recomputed the squared Frobenius norm from all 65,536 trial integers, verifying the bound 32. The new input-verification session executed the source oracle and kernel producer, checked the old positive pivots, and compared the exact targets without rerunning the sign. Both use the same implementations; separate execution is not independent implementation.

## Objects and comparison

All three matrix centers belong to the same actual even Legendre embedding

$$
E_j(u)=\sqrt{\frac{4j+1}{2a}}P_{2j}(u/a),\qquad j=0,\ldots,255,
\qquad a=\log3.
$$

The prime input supplies $A_0$ approximating $E^*M_9E$ and $B_0$ approximating **$E^*M_9^2E$**. The clipped prime-power shifts are $2,3,4,5,7,8$; the endpoint shift 9 has zero overlap. The second matrix is not $A_0^2$. The kernel center $J_0$ includes the polynomial kernel and the actual polynomial-filtered even rank-one tail update.

Every matrix center is an exact symmetric rational matrix on the $2^{-512}$ grid. The prime and regenerated kernel arithmetic operator allowances are at most $2^{-504}$. The kernel's analytic replacement uses the newly computed moment enclosures and the source's analytic proof. The consumer uses $\epsilon_J=2^{-205}$ to pay both that replacement and the directed finite assembly error.

For the exact rational trial $X$ on the $2^{-256}$ grid, its integer sum of squares verifies $\|X\|\le\|X\|_F\le32$. With $\mu=4/5$, the source inverse-residual supplier and matrix-error allowances give

$$
W=X+X^*-X^*A_0X+
\mu^{-1}(I-A_0X-X^*A_0+X^*B_0X)+\Delta_XI,
$$

$$
\Delta_X=1104\epsilon_A+1280\epsilon_B,
\qquad W\succeq G=E^*M_9^{-1}E\succ0.
$$

Thus $W^{-1}\preceq G^{-1}$, the direction required for a sufficient lower comparison. With $D_0=B_0-A_0^2$ and

$$
\epsilon_D=\epsilon_B+(2(2011/325)+\epsilon_A)\epsilon_A,
$$

the matrix tested by the consumer is

$$
T=W^{-1}+J_0-\epsilon_JI-2^{-187}I
-2^{-186}(D_0+\epsilon_DI).
$$

The common-embedding projection and complementary-block allowances $e=2^{-188}$, $n=2^{-390}$ and $h=2^{-386}$, and the coefficient comparisons for this target, are in the Library note's **A common 256-mode consumer** section. They are paper-level source applications, not kernel-verified declarations supplied by this package.

NumPy only proposes a possible negative vector; directed arithmetic checks its sign against $T$. Directed LDL supplies the recorded positive sign. The minimum recorded pivot lower bound is about $1.91066\times10^{-6}$; it is **not a minimum eigenvalue bound**. This fixed-width comparison does not supply the cofinal support family required for RH.

## Source inputs and regeneration

`input-manifest.json` binds the prime, kernel and moment bytes, scientific programs and selected environment. The `.b64` files contain base64-encoded zlib streams of their exact JSON files. Hashes identify bytes, rather than independently proving the execution or analytic claims.

The moment wrapper directly calls the pinned author's `complete_oracle()`, `check_rows()` and `mutation_controls()`. These functions require no author moment packet, retained matrix or certificate, and the moments do not depend on support half-width. Their definition is

$$
\vartheta_q=\frac{256}{\pi}\int_0^1x^{2q}
\bigl(\Re\psi(1/4+128ix)-\log\pi-7/2\bigr)\,dx,
\qquad 0\le q<1024.
$$

The complete Binet remainder is included once by the source oracle. Outward floor/ceiling rounding gives project intervals of maximum width one grid unit in the recorded executions, with the conservative midpoint error bound $2^{-1024}<2^{-1023}$ used by the analytic interface. The new rows replace the former author-packet containment premise; they do not certify the original author packet.

To regenerate just these moments and the affected kernel:

```sh
build/weil-c9-venv/bin/python docs/reports/weil-c9/c9_source_moments.py --checker build/liu-01-check.py --out build/fresh-c9-moments.json
build/weil-c9-venv/bin/python docs/reports/weil-c9/c9_kernel_ball_matrices.py --moments build/fresh-c9-moments.json --size 256 --bits 3072 --out build/fresh-c9-kernel.json
```

The prime input was computed with `c9_prime_ball_matrices.py --size 256 --bits 4096`; its recorded original runtime is 1795.43 seconds. It is reused without regeneration. For a separate inspection of that production:

```sh
build/weil-c9-venv/bin/python docs/reports/weil-c9/c9_prime_ball_matrices.py --size 256 --bits 4096 --out build/fresh-c9-prime.json
```

Runtime metadata can change raw output hashes while scientific fields agree. Reuse of prior sign evidence requires exact identification of all mathematical target inputs, errors and construction, as recorded here. Changed targets require their own sign evidence. An old execution must retain its original manifest and cannot be relabeled as a new execution.
