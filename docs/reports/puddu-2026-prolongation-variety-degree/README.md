# Prolongation-variety degree: exact checks

Source: Mattia Puddu, arXiv:2611.06638v1, Conjecture 4.1, and the first three jet pairs in Section 4. The source was submitted on 4 November 2026.

Run `python3 -u docs/reports/puddu-2026-prolongation-variety-degree/check.py` from the repository root. Dependencies: Python 3, SymPy 1.14, and Singular (tested with 4.4.1). The recorded output is `check.out`; its final line is `ALL_CHECKS_OK`. The program prints its own SHA-256.

The program builds the jets independently by the integer monomial rule for the derivation. It reproduces the source's first three jet formulas, the quadratic degree formula, and the cubic gcd and degree formulas on the tested range. It checks every exponent `2 <= n <= 5` and every order `0 <= m <= 12`, including both residue boundary cases. Rational root-isolating intervals verify simple negative roots and strict alternation of the reciprocal primitive factors through jet order 13. No floating point is used.

For every interior residue, both infinity charts are checked by SymPy Gröbner bases over the rationals. The equations use the actual homogenized jets and explicit integer coefficients on every lower jet. Adjoining `x^L,w^L` for `L=k+H+1` and `L=k+H+2` gives the basis `x^k,w^H` in both cases, where `k=n-r` and `H=(n-1)r`. There are 104 local Gröbner computations. These are independent computations on the equations, rather than counts taken from an assumed monomial ideal.

There are also 52 global Gröbner computations in Singular over the characteristic-zero coefficient field QQ, using degree reverse lexicographic order. At an interior residue the affine equations are

    f = A[n,m+1] + A[n,m+1-r]
    g = B[n,m+1] + 2 B[n,m+1-r].

At a boundary residue they are `f=A[n,m+1]`, `g=B[n,m+1]`. The program constructs these equations, runs `std(ideal(f,g))`, and obtains the quotient dimension with `vdim`. Every dimension agrees with the predicted degree, independently of the local-length subtraction.

These specialized sections can be nonreduced, so quotient dimension counts intersection length, not distinct points. The interpretation as the degree of a generic section additionally uses the mathematical argument: the two infinity unit conditions hold on a nonempty parameter open set, their local lengths are fixed there, and generic affine sections of the smooth embedded graph are finite and reduced. The finite computations do not establish the all-parameter conjecture. No Lean proof is included.
