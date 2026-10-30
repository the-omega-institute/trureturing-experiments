# Spezia A330700: higher exterior traces

`check.py` checks the second conjecture in the COMMENTS of
[OEIS A330700](https://oeis.org/A330700), revision 23, using the literal matrix
definition and examples in [A330287](https://oeis.org/A330287), revision 31, and
the perimeter interpretation in [A319840](https://oeis.org/A319840), revision 105.
The source indices start at one. Quotations and sequence values are from OEIS,
whose content is subject to the [OEIS End-User License](https://oeis.org/LICENSE).

Run `python3 check.py` with Python 3 and SymPy installed. The stored stdout is
`output.txt`; success means exit 0 and the final line `ALL_CHECKS_OK`.

The program reproduces the displayed matrices, eight published permanent
values, and geometric perimeter counts. It checks the proposed factorization
into a difference of two outer products through matrix order 40. Independently,
it computes exact dense characteristic polynomials through order 40, every
higher principal minor through order 10, and every entry of every higher
compound matrix through order 7. Determinants of restricted matrices use a
separate rational Gaussian elimination implementation. It also checks the
nonzero principal minor on source indices {1,3} and the zero-dimensional
exterior-space cases above the matrix order.

All mathematical computations use integers, exact SymPy algebra, or
`fractions.Fraction`. These finite checks support source fidelity and the
proposed settlement; they are neither a uniform proof nor Lean kernel
verification. The separate first conjecture's indexing mismatch is printed as
a scope diagnostic and is not the nominated statement.
