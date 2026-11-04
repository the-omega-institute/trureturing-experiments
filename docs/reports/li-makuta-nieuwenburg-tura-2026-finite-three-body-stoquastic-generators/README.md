# Finite three-body stoquasticity cone checks

Source: Jan Li, Owidiusz Makuta, Evert van Nieuwenburg and Jordi Tura,
*An exact semidefinite characterization of the stoquasticity cone of
permutationally invariant Bell operators*, arXiv:2610.03438v1, Definition 3,
Eq. (4), Section VI, and Section VIII's first open direction.
Predecessor: arXiv:2603.22493v1, Section V.2 and Section VII.

Run `python3 check.py` with Python 3 and SymPy 1.14.0. The saved output is
`output.txt`; success ends in `ALL_CHECKS_OK` and exit status 0.

All evidence uses integer or rational arithmetic and symbolic identities.
The seed 261003438 samples two pairs of rational half-angle tangents.
Their sine and cosine values are rational; these are exact angle settings.
For each pair and n = 3, 4, 5, 6, the program constructs the actual finite
Dicke-sector sums from ordered distinct-site tensor products, removes
lineality by pivot columns, enumerates all feasible one-dimensional active
faces, and compares the resulting complete ray set with the proposed one.
The computed quotient has rank 6 and respectively 6, 7, 8, 9 rays;
lineality has dimension 3.

Fidelity checks reproduce the source's two-body first- and second-off-diagonal
formulas and its three-ray, two-line description on these same cases.
Independent symbolic checks verify the proposed invertible off-diagonal
coordinate map and the collinear endpoint identity. Finite quadratic-cone
ray sets are also exhaustively checked for m = 2 through 12.

Collinear cases use four exact unit-circle settings, both signs of the
coincident measurement, and n = 3 through 6. They have rank 3, lineality
dimension 6 and exactly the proposed three quotient rays. Diagonal settings
have zero off-diagonal constraint matrices and the whole coefficient space
as lineality.

The unnormalized sector sums differ from normalized Dicke entries by positive
row-dependent factors. Consequently they define exactly the same inequalities.
All active-face enumeration takes place in the original Bell coefficients;
the proposed Pauli-coordinate reduction is checked separately against it.
These finite checks and symbolic identities support a preregistration draft;
they do not prove the universal generator theorem or certify a Lean proof.
