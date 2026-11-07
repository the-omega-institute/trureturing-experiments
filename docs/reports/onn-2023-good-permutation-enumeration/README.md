# Onn good-permutation enumeration

The target is the unnumbered referee question before Lemma 2.2 in Shmuel Onn,
*On Supmodular Matrices*, arXiv:2309.07700v1: count and characterize the rank
arrays whose adjacent squares have their maximum on the main diagonal and their
minimum on the opposite diagonal.

Run `python3 check.py`. The program uses Python 3's standard library and integer
arithmetic. Its output is `check.out`, ending in `ALL_CHECKS_OK`.

The fidelity checks reproduce the good rank arrangements in Theorem 1.2 and the
additional displayed examples. They also check all 55 nondecreasing nine-tuples
with values in {0,1,2} against the full supmodularity inequalities for the source's
3-by-3 arrangement. The source's Theorem 1.3 already excludes good arrays when
min(m,n) >= 3 and max(m,n) >= 4.

The independent enumeration covers every positive ordered dimension pair with
m*n <= 12. For two dimensions exceeding one, it explores the complete
permutation tree and prunes only when a completed square violates the literal
source condition. Each rejected prefix accounts for all factorially many
completions, and accepted plus rejected completions must equal (m*n)!.
For a single row or column, memoized exhaustive search over unused-label subsets
sums every permutation; it does not iterate separately through identical suffix
subproblems. Thus the 1-by-12 cases do not materialize 12! arrays.

The two-chain construction independently generates the candidate arrays. Exact
set equality with the enumerated good arrays checks both directions and
uniqueness of the split for widths 1 through 6, in both orientations. At width 6,
the split counts are 1, 66, 495, 924, 495, 66, 1, totaling 2048.
A separate unpruned loop checks all 9! three-by-three arrays against the two
proposed structural types: 36 of each, all with center 5. The 3-by-4 and 4-by-3
counts are zero.

These finite checks corroborate the proposed general argument. They are neither
a proof premise nor a Lean/kernel verification of the universal answer. The
journal version has not been compared in full; this report makes no novelty
claim.
