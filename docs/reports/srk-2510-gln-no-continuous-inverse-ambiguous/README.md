# Exact checks for the complex general linear candidate

Run `python3 check.py`. The program uses only the Python standard library, integer
arithmetic, and `fractions.Fraction`. Its output is stored in `check.out`.

Source: David Schmitz, Sadman Rahman and Anthony Kindness, *Continuous Inverse
Ambiguous Functions on Lie Groups*, arXiv:2510.09958v1, equations (1.1)–(1.2),
Theorem 1.3, section 2's even-torus construction, and the final paragraph of section
5 after Theorem 5.8; DOI 10.1007/s00010-024-01131-8.

The fidelity checks enumerate every bijection of the cyclic groups of orders
1 through 8 and verify Theorem 1.3's finite criterion and the equivalence of the
two defining equations there. They also check the even-torus formula on all
torsion points of exponents 2 through 9 in dimensions 2 and 4. These finite
restrictions test the literal formula, not its continuity on the full torus.

The settlement checks verify inversion/determinant identities on rational
invertible matrices in dimensions 2 through 8 and certify that no integer square
is minus one using all four residues modulo four. They check the rank-one
endomorphism composition law on small integers. No bounded sample supplies a
proof of a topological claim.

The proposed negative answer for every complex dimension greater than one still
requires a proved natural rank-one invariant of the actual matrix group. In
particular, a split surjection to the integers is insufficient: the program checks
an automorphism of the integer lattice of rank two whose square is minus the
identity and whose first-coordinate projection does not descend through that
automorphism. Neither the topology of the matrix group nor naturality for
arbitrary continuous self-maps is certified here.

The program is independently written from these mathematical formulas and
contains no imported third-party program or Lean source.
