# Equal-first-trace MM parent certificate

The mathematical statement and complete proof belong to the main
the-omega-institute/trureturing repository, relative path
docs/develop/theory/FIB_SELF_CALIBRATING_RELATION_MATHEMATICS_X.md, §14.
[result.md](result.md) describes the coefficient evidence and its scope.

Requirements: macOS or Linux, Python >=3.12, python-flint 0.9.0 and
SymPy 1.14.0 already installed. No repository checkout, network access or
installation is needed by these programs. Paths resolve from each script's
<code>__file__</code>; child interpreters use the invoking Python executable.

From this directory:

~~~sh
python3 -B run.py certificate
~~~

From any unrelated directory, including after copying this entire folder:

~~~sh
python3 -B /path/to/copied/fib_equal_trace_parent_180/run.py certificate
~~~

The entry constructs both universal **equal-first-trace** resultants with
FLINT, then independently reconstructs them with SymPy from the original
source quadrics. It checks exact divisions, the source determinant,
gcd(K,alpha)=1, the complete factor product, the coefficient-file hash,
and the degree-preserving residue 87 modulo 103. Each child has a
60-second deadline; the parent terminates and joins its exact child on
timeout and returns 124. Arithmetic errors return nonzero.

[same_trace_parent_F.json](same_trace_parent_F.json) is the primitive integer
polynomial in variable order (theta,xi); its bytes are preserved by verification.
[same_trace_verification.json](same_trace_verification.json) records the
exact-check results and resource measurements. Other generated JSON files
are mathematical inputs, complete coefficients or compact results.
[payload_manifest.json](payload_manifest.json) fingerprints every retained
file except itself; generated timing fields change when rerun.

Optional rational fixed-parent comparison:

~~~sh
python3 -B run.py fixed-images
~~~

It reproduces the two image curves and candidate supports at
(theta,xi)=(3/7,-35/29). Their common support is 23b-29 up to a nonzero
constant. This comparison is not used to prove the universal theorem.
The retained model also exposes the entry run.py universal-trial, a bounded
search for the remaining unequal-first-trace incidence. It supplies no global
certificate and is not part of certificate reproduction.

Exact arithmetic, the prose source/world and continuation proof, independent
mathematical review, Lean verification, CI and publication are different
evidence categories. The scripts establish only the arithmetic category.
