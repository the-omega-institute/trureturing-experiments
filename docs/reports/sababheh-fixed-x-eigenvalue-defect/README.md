# Fixed-X higher-eigenvalue defect: exact diagnostics

Source: Mohammad Sababheh, *An Extremal Formula for Elliptical Width and the
Bourin–Lee Conjecture*, arXiv:2610.04377v1, Lemma 2.1, Theorem 3.2, and the
unnumbered fixed-X question in Section 6. Eigenvalues are decreasing and one-based;
elliptical width uses the **diameter**, rather than radius, of an inscribed disk.

Run with Python 3 and SymPy 1.14.0:

```sh
python3 check.py > check.out
```

The program uses rational and Gaussian-rational arithmetic throughout. Hermitian
inertia is computed by exact congruence with scalar or 2-by-2 Schur pivots; spectral
rank bounds count eigenvalues above rational thresholds. It does not use floating
point, approximate roots, or numerical eigenvalue routines.

The source-fidelity checks reproduce the disk support function in Lemma 2.1 for
`beta I + a E12` symbolically, and Theorem 3.2's norm-defect value for three
nilpotent examples with explicit maximizing positive completions. The supporting
characteristic polynomial has roots `Re(u beta) +/- a/2` for unit complex `u`;
convexity of the numerical range therefore gives the disk of radius `a/2`.
The first two coordinates attain that disk, so its elliptical width is `a`.

For every dimension 2 through 6 and every allowed higher index, the program checks
90 positive definite completions for the invertible nonnormal matrices
`X = (3+4i) I + a E12`, where `a` is 1 or 2. At interior indices their defects lie
between `a-error` and `a`, with `error < 1/q`; terminal defects lie between
`-eps-eta` and zero, with `eps+eta < 2/q`. The exact budgets are `q = 2,4,8`.
An additional 15 cases retain nonzero cross-compression blocks and check positivity
and both spectral ranks in the proposed general lifting lemma.

These are finite exact diagnostics of a general proposed argument. They do not
compute the unrestricted supremum by enumeration and do not establish its value
for arbitrary matrices or dimensions. The proposed complete answer is elliptical
width for `1 <= j <= n-2` and zero for `j = n-1`. Its proof still requires the
general compression-lifting argument and the published upper bounds.

Producer: independent codex-cli research seat, flight `424be8766a026ed9e12c516f`;
prior disclosure `repo-prior-exposed`. No Lean proof is part of this report.
