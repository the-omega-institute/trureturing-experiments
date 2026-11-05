# Cazenave Conjecture 5: exact finite check

Source: Tristan Cazenave, *Lemmas and Conjectures in Spectral Graph Theory*, arXiv:2611.03879v1, Section 7, Conjecture 5 (corpus 239534). The statement is μ₂(L) + rank_R(A) ≤ 2m for every connected simple graph of order n ≥ 2. Eigenvalues of the ordinary Laplacian are ordered decreasingly, with multiplicity.

Run with Python 3, SymPy, and nauty 2.8.9:

```sh
python3 check.py --geng /path/to/nauty2_8_9/geng
```

Generator source: https://pallini.di.uniroma1.it/nauty2_8_9.tar.gz, SHA256 `c97ab42bf48796a86a598bce3e9269047ca2b32c14fc23e07208a244fe52c4ee`. Extract it and build with `LC_ALL=C ./configure` followed by `LC_ALL=C make geng`. The check invokes `geng -cq n` for every n from 2 through 9 and independently verifies connectivity, graph6 decoding, and the known connected unlabelled graph counts (1, 2, 6, 21, 112, 853, 11117, 261080). The exhaustive scope is up to isomorphism; the target is invariant under vertex relabelling.

All arithmetic used as evidence is exact. Integer symmetric congruence computes the adjacency rank r and the inertia of (2m−r)I−L. Its negative inertia counts the Laplacian eigenvalues strictly above 2m−r; at most one is equivalent to the desired bound on the second largest eigenvalue. A nonzero diagonal gives a one-dimensional Schur complement. When all diagonals vanish, a nonzero off-diagonal gives a two-dimensional pivot with one positive and one negative sign. Every remainder is scaled by a positive integer and divided by its positive common gcd, preserving inertia. This method handles singular matrices and repeated eigenvalues without numerical tolerance.

The program cross-checks this routine against SymPy's characteristic polynomials and rational root isolation for all 729 symmetric 3-by-3 matrices with entries in {−1,0,1}, and against exact ranks and root isolation for every generated graph with n ≤ 5. Fidelity checks reproduce the source's equalities for K₂, K₃ and P₄ and the spectra/ranks of the other trees of order at most four.

The stored output covers all 273192 connected unlabelled graphs with 2 ≤ n ≤ 9, with no failures. The exact minimum slack 2m−rank_R(A)−μ₂ is 0, attained by K₂, K₃ and P₄. No other zero slack is observed in this finite range. This is a finite computational check, not a Lean proof, a claim about graphs with n > 9, or a global equality classification. The proposed positive settlement uses a general centroid and spectral argument; this enumeration is supporting evidence.
