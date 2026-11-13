# Exact checks for the higher-rank constant-width conjecture

The target is the unnumbered Section 1 conjecture of Chien–Li–Nakazato,
*The diameter and width of higher rank numerical ranges*, Linear and Multilinear
Algebra 69(5), 871–887, DOI 10.1080/03081087.2019.1710105.
The source manuscript is <https://cklixx.people.wm.edu/diameter.pdf>.
Jiang, arXiv:2611.13461v1, treats the classical case and explicitly excludes
this higher-rank extension.

Run `python3 check.py` with Python 3.12 or newer and SymPy 1.14.0.
The stored output is `check.out`. Every comparison uses rational arithmetic in
Q(√3), or exact SymPy algebra; floating point is not used as evidence.

The source fidelity check reproduces the symbolic Hermitian-block characteristic
polynomial and tests Theorem 3.1's eigenvalue gap at 184 rational unit directions
and four admissible radii. The proposed counterexample uses the source's own
A(√3), with rank 3 in dimension 6. Its proposed higher-rank range is the
intersection of the radius-√3 disks centered at 1, ω, ω², a Reuleaux triangle.

For each tested direction the program exhibits an attaining boundary point and
a global support certificate: containment in a defining disk, or a nonnegative
combination of the two active disk normals at a vertex. All 184 exact widths
are √3. The interior-ball inequality and the obstruction to being any disk are
exact. Another 4,416 cases check the eigenvalue threshold for ranks 3 through 8
with zero through three centroid scalar blocks. Eight explicit algebraic
orthogonal projections verify scalar compressions at the centroid and the three
vertices, at ranks/dimensions 3/6 and 4/9, including affine transport.

These are finite exact checks and symbolic identities. They do not replace the
universal disk-intersection and constant-width arguments in the preregistration,
or claim a Lean proof. The intended projection definition is PAP = zP; the
manuscript's AP = zP conflicts with its own rank-1 assertion and its quoted
Li–Sze characterization. No Lean claim or negation is proved here.

The general family still requires a formal proof. The cases k = 2 and n < 2k
are outside the proposed family. The classical k = 1 conjecture is a different
statement. Source verification and literature screening belong to the draft;
no originality claim is established by this computation.
