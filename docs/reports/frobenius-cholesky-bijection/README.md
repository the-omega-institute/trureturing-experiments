# Frobenius-compatible triangular zero fibers

Source: Hongfeng Wu and Kai Zhou, *Singular Cholesky Fibers over Finite Fields*,
arXiv:2609.00533v1, Sections 6.1 and 8(iv).
https://arxiv.org/abs/2609.00533v1

Run with Python 3 and its standard library:

    python3 check.py > output.txt

The fidelity check reproduces the even-characteristic rank recurrence for
square-zero triangular matrices and Cholesky roots of zero. It independently
enumerates both families for every size 0 through 5 over GF(2), 0 through 4 over
GF(4), and 0 through 3 over GF(8). Their rank distributions agree with
q^r f(n-1,r) + (q^(n-r)-q^(r-1)) f(n-1,r-1), with f(n,0)=1.

For every square-zero input in these ranges, the program checks the proposed
minor-based recursive map, its inverse, upper triangularity, the zero Gram
condition, rank preservation, Frobenius compatibility and leading-principal
truncation. Surjectivity and both inverse identities are checked against an
independent enumeration of the entire Cholesky-root family. Both embeddings
GF(4) into GF(16) are checked on all 496 size-4 square-zero inputs. The last
output line is ALL_CHECKS_OK.

Arithmetic is exact in binary polynomial quotient fields with moduli
x+1, x^2+x+1, x^3+x+1 and x^4+x+1. Addition is integer XOR; multiplication and
inversion use exact tables. Field laws are checked exhaustively. Determinants,
minor selection and Cramer's rule implement the construction. A separate
Gaussian-elimination rank checker validates the resulting ranks. Enumeration
of integer field labels belongs to the test harness; the map's choices use
only coordinate indices and nonzero-minor tests.

These finite checks do not prove the assertion for unbounded size or all
characteristic-two fields. The planned positive settlement is a general
inductive argument, not bounded enumeration. This report contains original
check code and its output, and does not redistribute the paper.
