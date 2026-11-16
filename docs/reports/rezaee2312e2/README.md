# Rezaee Conjecture E(2): exact checks and source convention boundary

Source: Fatemeh Rezaee, *Conjectural criteria for the most singular points of
the Hilbert schemes of points*, arXiv:2312.04520v2, Section 1.4, Conjecture E(2).
The arXiv metadata identifies v2 as the revision of 2026-09-25. The Cambridge
repository copy is labelled “Accepted version”; the publisher's version of
record has not been read.

Run `python3 docs/reports/rezaee2312e2/check.py` with Python 3.12 or newer and
Singular on PATH. The stored output was obtained with Python 3.14.8 and
Singular 4.4.1. All mathematical computations use integers or characteristic
zero exact algebra. The final line `ALL_CHECKS_OK` means the stated checks
passed; it does not certify the conjecture's global maximality or geometry.

The script solves the full coefficient system for pairwise monomial
syzygies by equality classes with a forced-zero class. It independently
computes the same Hom dimensions using Singular's syzygy module over
the Artinian quotient. Scalar extension from rational numbers to complex
numbers preserves the ranks of these rational matrices.

The fidelity checks reproduce the source's formula for the powers of the
maximal ideal and its explicit E(2) example of colength 12 and tangent
dimension 66. Both algorithms give these E(2) dimensions:

| k | Colength | Tangent dimension of I_k and J_k | Borel staircases with smallest pure exponent k |
|---|---|---|---|
| 1 | 3 | 9 | 2 |
| 2 | 6 | 24 | 2 |
| 3 | 12 | 66 | 2 |
| 4 | 22 | 156 | 2 |
| 5 | 37 | 321 | 2 |
| 6 | 58 | 594 | 2 |

Here J_k has standard monomials of degree less than k, together with z^k
and yz^(k-1). Its homogeneous tangent dimensions match the proposed
degree 0, -1 and -2 formulas in every checked case. The Borel enumeration
uses the two additional standard monomials above degree k-1 and verifies
order-ideal and strong-stability-complement conditions.

The source notation defines conv(I) using the exponents of **all monomials
of I**, whereas Conjecture A describes a hull spanned by its monomial
generators and prescribes an upper simplex. These hulls differ already on
every checked I_k: z^(k+3) belongs to I_k but lies outside the compact hull
of the minimal generators, whose z-coordinates are at most k+2. Both
interpretations cannot be treated as one literal definition. This diagnostic
requires clarification before encoding the entire geometric assertion; it
does not refute the intended compact-generator-hull conjecture.

No global maximality, flat-family theorem, universal syzygy count, or convex
geometry clause has been formally proved here. No Lean proof is included.
