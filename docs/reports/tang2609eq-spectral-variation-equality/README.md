# Equality in Tang's truncated spectral-variation bound

This program checks the complete proposed answer to Qiyue Tang, *A Truncated
Singular-Value Bound for Spectral Variation of Normal Matrices*,
arXiv:2609.09177v2, Section 3, Question 5. The source inequality is

\[
d(\sigma(A),\sigma(B))^2\leq
\sum_{j=1}^{\lfloor(n+1)/2\rfloor}s_j(A-B)^2.
\]

The proposed characterization has three branches: equality always in dimension
one; simultaneous diagonalization with an optimal direct pairing in dimension
two; and, in dimensions at least three, either identical matrices or one
different common reducing eigenline, with no path of shorter displacements
connecting its two distinct eigenvalues through the common spectrum.

Run `python3 check.py` with Python 3 and SymPy 1.14.0. The stored result is
`output.txt`. The program uses integer and rational arithmetic, exact complex
conjugation, characteristic polynomials, and certified rational intervals
isolating their real roots. Strict inequalities are accepted only when rational
interval endpoints separate them. Equality is accepted only when the relevant
interval is a singleton. No floating-point computation supports a conclusion.

The fidelity check reproduces the source example
A = diag(t,0,…,0), B = 0, at t = 1/3 and t = 2 in dimensions one through four.
The program tests every pair of diagonal matrices over {0,1,2} in those
dimensions. It also tests rational-entry normal matrices in dimensions two,
three and four: independent rational unitary bases with Gaussian-integer
spectra, common bases, repeated eigenvalues, noncommuting rank-one differences,
rank-two equality in dimension two, strictness after rematching, and strict
pairs close to equality. Random choices use the fixed integer seed 260909177.

Matching distances are computed by exhaustive permutations of the spectral
lists. The structural predicate is computed separately: joint eigenspace
restrictions in dimension two, and the range of the difference plus a labelled
threshold graph in higher dimensions. The full finite tests support the draft;
they do not prove the assertion for arbitrary complex normal matrices or any
unbounded dimension. The proposed Lean settlement uses a general argument;
these finite computations are auxiliary checks.
