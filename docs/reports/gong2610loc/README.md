# Exact pinned-residue locality checks

The source is Gong, Lyu, Hu, Li and An, *Correlations decide a shallow-circuit
advantage*, arXiv:2610.04666v1, Section V(c), Eq. (6). The proposed result removes
the additional fanout and shared-seed restrictions for exact residue pinning:

N <= 4^max(d-2,0) sum_{j=1}^{min(k,d)} binom(k,j).

Run `python3 check.py`. The standard-library program uses integers and
`fractions.Fraction` throughout. Its complete output is in `output.txt`.
No third-party program or Lean source was copied.

- Fidelity: with three nonempty Walsh characters on two read seeds, one unused
  seed and p=3, each of the three signed fiber vectors gives Eq. (6)'s weighted
  sum -4, congruent to -2^(3-1) modulo 3.
- Exhaustive scope: all Boolean truth tables on ell=0,1,2,3 fair seed bits; all
  unordered balanced, pairwise-orthogonal output families, including all output
  polarities; all read subsets K; d=0,1,2,3 whenever the outputs are d-local;
  p=3,5,7,11. The 63,162 families yield 262,544 parameter cases. Every nonempty
  family has N<=7 by enumeration. Prime 11 also represents every odd prime
  greater than 7, since equality of weights modulo such a prime is equality
  of weights as integers on this finite range. Output ordering changes none
  of the tested statements, and repeated or complementary outputs cannot
  occur in an orthogonal family.
- Exact ingredients: 1,226,536 source congruence checks, 505,163 expanded
  sum-of-squares/Bessel checks, Fourier support and granularity checks, and
  the projected-energy lower bound at each pinned configuration.
- Sharpness: all nonempty Walsh characters of degree at most d on K, for
  d=1,2 and k=0,...,6, attain the first bound. Pinning is checked at
  p=3,5,7,11. The general sharpness argument is in the preregistration draft.
- The parameter calculation d=k=3 gives N<=28.

The final output is `ALL_CHECKS_OK`. These finite checks do not establish the
unbounded theorem; the proposed settlement uses a general divisibility,
Fourier-granularity and Bessel argument. No Lean proof of the target or its
negation is included.
