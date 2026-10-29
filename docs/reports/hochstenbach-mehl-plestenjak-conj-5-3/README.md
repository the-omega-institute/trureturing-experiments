# Hochstenbach–Mehl–Plestenjak Conjecture 5.3: exact checks

Source: *Solving singular generalized eigenvalue problems. Part III: structure preservation*, arXiv:2406.07109v1, Section 5. The claim concerns real symmetric singular pencils and real symmetric rank-completing updates of arbitrary rank.

Run `python3 check.py` with Python 3 and SymPy. The accompanying `output.txt` was produced with Python 3.14.8 and SymPy 1.14.0. All matrices, determinants, discriminants, resultants and root-order checks use exact arithmetic. Successful execution exits 0 and ends with `ALL_CHECKS_OK`.

The program reproduces the normal rank two and minimal index one of the source's displayed 3×3 pencil, its prescribed eigenvalue 2, and the double random eigenvalue with geometric multiplicity one and the stated orthogonality. There are two sign discrepancies in the source's numerical examples: the displayed 3×3 pencil gives random eigenvalue −1 rather than +1; the displayed 5×5 pencil gives the factor λ²+λ+1 rather than λ²−λ+1. The program records both explicitly. These discrepancies do not contradict Conjecture 5.3.

The proposed determinant identity is checked as a polynomial in both λ and τ for every ordered minimal-index list of length 1–3 with entries 0–2 and regular-part size 0–2, subject to total pencil size at most 7: 42 size choices, with three specified update matrices each, giving 126 identities. The prescribed pencil is symmetric and includes off-diagonal entries. An additional size-four case keeps all eight update entries symbolic. These finite checks do not prove the identity for arbitrary sizes or matrices.

A size-nine, rank-two witness has observation determinant (λ−4)(λ−5)(λ−6), discriminant 4 and nonzero resultant with the fixed spectrum. Further checks cover two eigenvalues at infinity (one true and one prescribed), repeated coincident fixed eigenvalues, zero minimal indices and an empty regular part. No Lean proof or kernel verification of the conjecture is supplied here.
