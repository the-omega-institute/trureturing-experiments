# Hanna A359914: exact partition-parity checks

Run `python3 check.py`. The standard-library-only program uses Python integers
throughout and exits zero with final line `ALL_CHECKS_OK`. Its output is in
`output.txt`.

The program solves the original bilateral equation by coefficient recursion,
using the two distributed source monomials at integer index `k`:

```
(-1)^k x^(k(3k-1)/2) A^(3k)
+ (-1)^k x^(k(3k+1)/2) A^(3k+1).
```

Both shifts are at least `abs(k)`, so the degree-`N` coefficient uses exactly
the finite window `-N <= k <= N`. Negative powers of `A` are computed with an
exact integral unit inverse. These are the same summand and coefficient window
as the preregistration's binding Lean conventions; the computation does not
replace the defining equation with the proposed congruence.

Results through degree 40:

- All 23 listed source terms and the b-file coefficients at indices 0 through
  40 agree exactly with the recursively generated coefficients.
- Every coefficient of the original equation, an enlarged window, the
  reindexed paired equation, and source FORMULA (2)'s Watson product agrees.
- An independent multiplicity dynamic program counts ordinary partitions.
  Its series is reciprocal to the checked Euler product.
- At every index 1 through 40, `a(n)` is even and `a(n)/2` has partition parity.
- The proposed modulo-four pairing also holds on 20 additional integral units
  congruent to one modulo two.

These are finite fidelity and consistency checks. The all-index argument is a
preregistration outline, and no Lean proof of the claim or its negation is
part of this experiment.

Sources: Paul D. Hanna, OEIS A359914, authored January 23, 2023, revision #43
(January 19, 2024, 10:06:08), first COMMENTS conjecture, NAME, FORMULA (1) and
(2), and EXAMPLE. Immutable record:
<https://raw.githubusercontent.com/oeis/oeisdata/571197984aeda7b8b680fe3a1f6a43f9898cd441/seq/A359/A359914.seq>.
Record SHA-256: `2706b4b0d1acf358b0cc1e2eb0a1a73378e12e4d3602fbe7d0972988c92cef89`.
The downloaded <https://oeis.org/A359914/b359914.txt> has SHA-256
`f65c64db4fe627b98d753f4269e1f3ea86fb19bd9b24ff5d2164ea95d5752a0a`.
The numerical prefix embedded in the checker is attributed to Hanna and the
OEIS Foundation; OEIS data are available under CC BY-SA 4.0.

The FORMULA introduction's `Sum_{n>=1}` omits the constant term. NAME,
FORMULA (1), the displayed EXAMPLE, and the numerical data require `A(0)=1`.
The checker uses that reading. The second COMMENTS comparison with A361552
is outside this experiment's selected assertion.
