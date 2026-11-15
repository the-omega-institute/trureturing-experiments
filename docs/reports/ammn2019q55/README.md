# Solomon–Terao divisibility converse

Exact arithmetic checks for Question 5.5 of Abe–Maeno–Murai–Numata,
*Solomon–Terao algebra of hyperplane arrangements*, JMSJ 71 (2019),
1027–1047, DOI 10.2969/jmsj/79957995, arXiv:1802.04056v1.

Run `python3 check.py` with Python 3.12 or newer and SymPy 1.14.0.
The output is in [check-output.txt](check-output.txt) and ends with
`ALL_CHECKS_OK`. The program uses integers and exact rational polynomial
algebra throughout.

The fidelity check reproduces the polynomial in source Example 5.4 for
`xyz(x+y+z)` and the displayed Boolean deletion polynomial. For the proposed
family with q = 4, 5, 7, 8 it checks distinct normals and the ordinary restriction,
logarithmic bases and their Saito determinants, admissibility of the sum of
squares on every intersection flat, the Hilbert polynomials from Gröbner bases
of the explicit Solomon–Terao ideals, Hilbert divisibility, and failure of
Poincaré divisibility from an independent Möbius lattice calculation.

Pure-power leading monomials bound every standard monomial, so the Hilbert
calculation exhausts all degrees. Finite samples support the stated infinite
family; the family proof is a separate mathematical argument. These checks
are not Lean proofs or a claim of kernel verification.

The check program is independently written; SymPy is used as a dependency
under its BSD licence. No program source from a settlement repository is used.
