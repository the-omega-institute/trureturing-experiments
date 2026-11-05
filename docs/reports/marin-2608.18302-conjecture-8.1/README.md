# Twisted Schur specializations with two dominant Newton vertices

The target is Carles Marín, *Schur polynomials twisted by roots of unity and
reciprocal pairs: torsion filters, fusion quotients, and total unimodularity at
odd order*, [arXiv:2608.18302v2](https://arxiv.org/abs/2608.18302v2), Conjecture 8.1.

Run `python3 check.py` with Python 3 and SymPy (verified with SymPy 1.14.0).
The stored standard output is `check-output.txt`; success ends in `ALL_CHECKS_OK`.
Every calculation uses integer, rational, or exact symbolic arithmetic.

The program reconstructs the elementary functions of
`(1, omega, omega^2, z, z^-1, w, w^-1)` modulo `omega^2 + omega + 1`.
It constructs the dual Jacobi–Trudi matrix for `(4,4,1,1)`, checks all 33 Laurent
coefficients, and computes the convex hull with integer cross products. Its
vertices are the signed permutations of `(4,0)` and `(3,3)`. Coordinate and
coordinate-sum maxima would force `(4,2)` for a single orbit; the integer normal
`(3,1)` separates that point from the hull.

As a fidelity check, the program independently expands the denominator
alternant and reproduces the source's denominator-polytope proposition at
`t=3, r=2`: its hull is `conv W(6,4)`. Frozen 3-by-3 Fourier minors are checked
for every residue triple; their common determinant has square `-27`.
Laplace expansion then computes the seven-variable alternants over the integers.
The base numerator equals the product of the 33-term polynomial and denominator.

For `b=(10,9,5,4,2,1,0)` and `s=1,4,7`, the program checks that the numerator at
`s*b` is the base numerator with both exponents multiplied by `s`. It verifies
three exact rational evaluations of the literal alternant quotient per `s`,
and checks the support subtraction on 90 chamber directions per `s`.
The exposed dominant vertices are `(10*s-6,4*s-4)` and `(9*s-6,7*s-4)`.
These finite checks support the family argument; they are not a Lean proof of
the conjecture's negation or a verification of every `s`.
