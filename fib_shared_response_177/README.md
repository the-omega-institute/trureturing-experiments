# Fixed-parent whole-factor separation

This bundle supplies the exact arithmetic for section 11 of
`docs/develop/theory/FIB_SELF_CALIBRATING_RELATION_MATHEMATICS_X.md` in
`the-omega-institute/trureturing`. It consumes the unchanged
`fib_shared_response_176/global_preflight` coefficient files at experiment
commit `7c45a21dbea5cdc0974f878a96dc64281cc5db07`. The mathematical result is
the elimination of three whole factors from the necessary support of the
**joint** independent MM/MJ collision question at the fixed parent
`(theta,xi)=(3/7,-35/29)`.

Python >=3.12, SymPy 1.14.0 and python-flint 0.9.0 are required, as pinned
by the consumed bundle's `requirements.txt`. With that environment, run:

```sh
python3 /path/to/experiment-repository/fib_shared_response_177/verify.py
```

The entrypoint resolves inputs relative to its own file. It accepts any
working directory and paths containing spaces; it uses neither network
access nor an author-specific absolute path. It reads the old coefficient
files without executing their producers or writing into their directories.
The default command recomputes and compares `certificate.json`; `--write`
regenerates that mathematical certificate beside the entrypoint.

The default output has 105 successful checks (104 mathematical/environment
conditions and one comparison with the preserved certificate), with
`MM_necessary_roots=104`, `joint_necessary_roots=54` and
`new_joint_exclusions=57`. Exit zero certifies these concrete computations;
it is not a repository gate or independent review.

| Block | Exact certificate and consequence |
| --- | --- |
| F44, old `support_factor_7.json` | Irreducible degree 44, eight distinct real roots. Its valuations in `Res_u(R,R_u)` (degree 124) and `Res_e(H,H_e)` (degree 240) are both one. It is coprime to the leading coefficient, image content and all A/D/Phi chart exceptional resultants. Hence R and H each have exactly one double geometric root, and the resultant image multiplicities force MM to separate all distinct R roots. |
| F0, `23b-29` | At `b=29/23`, R is reduced of degree 16 and `gcd(R,N)=gcd(R,N*PhiSharp)=u^2-u-1`; PhiSharp is a unit. Those two roots have `n=N/D=0`, so neither can be a positive world. On the remaining reduced degree-14 source, the MJ denominator is a unit; the exact cyclic determinant is nonzero, certified by residue 43 modulo 101. |
| F91, old `support_factor_8.json` | Irreducible degree 91, 49 distinct real roots; R is reduced and all chart exceptions are units over `Q[b]/F91`. The exact homomorphism `b -> 41` to `F_101` gives a reduced degree-16 source, an invertible MJ denominator and cyclic determinant 47. This certifies nonvanishing in the number field, and therefore separation in every real embedding. |

The finite-field statements are deterministic nonzero certificates, not
probabilistic modular tests or samples of real worlds. All coefficient
denominators and the source leading coefficient are invertible at the
displayed homomorphism. The supplied Bezout identities invert the MJ
denominator; consequently the source algebra and its cyclic matrix reduce
at that homomorphism. A zero characteristic-zero determinant would reduce
to zero, contradicting 43 or 47. For F91, irreducibility then carries the
nonzero algebraic number to a nonzero number in every embedding. The
Vandermonde identity proves that a full-rank matrix with columns
`1,q,...,q^(m-1)` is equivalent to separation of the distinct source roots.
The program checks the concrete cyclic determinant with both FLINT and an
independent integer Gaussian calculation modulo 101.

The program also reconstructs the VI24 adjugate and Sylvester identities,
the chart and MM image identities, actual cumulative M/J traces, the entire
cubic interpolation and its parent/sixth/seventh-clock identities, the
actual MJ numerator/denominator identity, all chart units, irreducibility
and exact Sturm counts. `certificate.json` holds input hashes, exact
modular polynomials, Bezout coefficients and cyclic matrix columns in
ascending powers of u. Its F0 polynomials use `b=29/23 -> 32 mod 101`;
its F91 polynomials use `b=41 mod 101`. The stored U and W need not be
reduced modulo the displayed source polynomial.

The MM necessary support loses only F44: degree 202 becomes 158, and 112
real roots become 104. The joint necessary support also loses F0 and F91:
degree 66, with the 54 distinct real roots of old factors 1--6. X9 already
excludes one F91 root in its Bprime band, verified here by exact endpoint
Sturm counts; it is not counted as a new result. There are 57 new joint
candidate exclusions relative to that existing 111-root baseline.

The complete-world correspondence is VI24.5: the real roots are filtered
by the original Q/T/child-line equations, strict entrywise positivity and
Phi != 0, then the unique drift and entire gain polynomial are recovered.
The finite complex source sets used here are supersets of these worlds;
no assertion that each algebraic root is a legal world is needed. No
Delta, sampled gain, highest gain coefficient or sign is divided out. The
two collision pairs remain independently quantified on the same actual b.
The 54 residual roots include all remaining chart/leading/content
exceptions and require full VI24 specialization; they are not asserted to
be actual collisions. Fixed-parent joint emptiness, the other sixth arm,
whole-parent and ancestor obligations, the original-root eight/nine-read
decision and a nine-read upper bound remain unresolved.

The ordinary algebra proofs are in X11. Resultants, multiplicity,
Vandermonde matrices, Gauss's lemma and Sturm counting are standard tools;
the fixed FIB block exclusions and their full-world application are
repository-derived. These certificates are not Lean proofs. Source code
and coefficient data are covered by the repository's Apache-2.0 license
in `../LICENSE`; no third-party source is copied into this bundle.
