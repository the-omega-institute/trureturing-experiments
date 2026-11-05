# Chapoton A277506: maximal chains of Hochschild lattices

Run `python3 check.py` with Python 3. No third-party packages are required.
All arithmetic is exact integers or `fractions.Fraction`.

The source is F. Chapoton’s conjecture in OEIS A277506, revision 24, and
Camille Combe, *A geometric and combinatorial exploration of Hochschild
lattices*, arXiv:2007.00048v1, Sections 1.1–1.2 and the longest-chain lemma
(journal version: Electron. J. Combin. 28(2), P2.38, DOI 10.37236/9929).

The program reproduces Combe’s explicit triword sets for dimensions 1–3,
the vertex formula and longest-chain length for dimensions 1–8. It computes
covers directly by testing all componentwise strict intervals for an
intermediate vertex. Dynamic programming over this independently derived
Hasse graph counts every inclusion-maximal chain. Dimensions 1–8 give
1, 2, 8, 51, 460, 5425, 79206, 1377985, matching the first eight OEIS terms.
The pentagon has chains of both lengths 2 and 3; both are counted.

The proposed cover classification and two-parameter continuation recurrence
are checked at every relevant vertex. Exact rational series calculations
through degree 12 check the defining Lambert W equation, the e.g.f., and the
OEIS recurrence. Stored output is `check.out`; its last line is
ALL_CHECKS_OK. These finite checks do not prove the conjecture for all sizes.
