# Hardin A126473 alphabet increments

For nonnegative alphabet size b and positive length n, A(b,n) counts words on
{0,...,b-1} whose adjacent absolute differences are at most 3. Leading zeroes
are allowed. The nominated empirical identities are
A(b,n) = A(b-1,n) + 7^(n-1) for b >= 3n-2, and
A(b,n) = A(b-1,n) + 7^(n-1) - 2 for n >= 2 and b = 3n-3.

Source: R. H. Hardin, [OEIS A126473](https://oeis.org/A126473/internal),
first COMMENTS line, NAME and S/R program; revision 57, 13 March 2026.
The 24 published sequence terms used for fidelity are attributed to OEIS;
[OEIS terms](https://oeis.org/LICENSE) apply to the source database.

Run with Python 3, without external packages:

    python3 docs/reports/hardin-a126473/check.py

The stored output is `check-output.txt`. Arithmetic is exact Python integer
arithmetic. The checker reproduces all 24 source terms at alphabet size 5,
compares dynamic programming with literal word enumeration in 45 cases,
checks both empirical clauses for n=1..12 and b=1..3n+4, and independently
checks translation fibers by exhaustive increment enumeration for n=1..7.
It also checks the proposed general step-bound formulas for radius 1..4 and
length 1..10. The final line is `ALL_CHECKS_OK`.

The checks are finite evidence. The proposed uniform proof uses the number
of translates of each increment word and classification of the two words
of maximum possible span. Neither this program nor its output proves the
unbounded assertion, and no Lean proof is included here.
