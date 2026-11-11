# Qi–Cui–Xu: longer dependency cycles

The proposed family has one edge at `(k+1,1)` and two-edges
`((i,1),(i+1,2))`, with row successor interpreted cyclically on `1,…,k`.
The associated real bilinear SOS has displayed length `k+1`.

Run `python3 docs/reports/2607.07071-section-9-problem-2-long-cycles/check.py`
with Python 3.12 or later and SymPy 1.14. The stored output is `check.out`.
All evidence uses exact integer, rational or real algebraic arithmetic.

The fidelity checks reproduce the source's Section 5 graph (12 edges, one
complementary dependency pair, five pairs subject to W3) and Section 8
Gram identity, exact spectrum and rank drop from 11 to 10. The family checks
cover every integer `3 ≤ k ≤ 24`: simplicity, C4-freeness, W2', W3, exact
dependency arcs, absence of complementary pairs, a directed k-cycle and
connectivity on the occupied row/column vertices. Specializing at `y=(1,0)`
gives the identity quadratic form of rank `k+1`.

For `3 ≤ k ≤ 6`, SymPy factorization over Q and Q(√2) gives one polynomial
factor. The real-field conclusion uses an additional exact certificate:
the quadratic coefficient matrix in x specializes to the identity, so its
rank over C(y1,y2) is at least three; the coefficients y1² and y1²+y2² have
gcd one. A quadratic product of linear factors has rank at most two, and
primitivity excludes a factor in y alone. Thus these checks also certify
irreducibility over R (indeed absolute irreducibility). Q(√2) is not claimed
to be the whole real field. Polynomial factorization is distinct from the
source's meaning of irreducibility, which is minimal real bilinear SOS length.

The all-k proposed argument is the specialization rank bound: every r-square
bilinear decomposition specializes to a linear-square representation of the
identity quadratic form, yielding MᵀM=I and r≥k+1. The finite checks support
this argument; they are not a Lean proof or an exhaustive positive settlement.

Source: Liqun Qi, Chunfeng Cui and Yi Xu, *A Weak Condition for Limited
Augmented Zarankiewicz Numbers*, arXiv:2607.07071v4; Mathematics 14(17),
3052, DOI 10.3390/math14173052. W3 and opposite-cell conventions follow the
corrected final journal definitions.
