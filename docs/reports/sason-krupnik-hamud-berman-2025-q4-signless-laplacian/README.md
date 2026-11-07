# Signless Laplacian spectra of the four subdivision–R partial joins

Source: I. Sason, N. Krupnik, S. Hamud and A. Berman, *On Spectral Graph
Determination*, Mathematics 13(4), 549, DOI 10.3390/math13040549, Question 4;
arXiv:2412.20775v6, Question 5.29. The underlying constructions are from
A. Das and P. Panigrahi, DOI 10.46793/KgJMat2306.947D.

Run `python3 check.py` with Python >= 3.12 and SymPy 1.14.0. The stored
`output.txt` ends with `ALL_CHECKS_OK`; the program exits 0. All matrices,
characteristic polynomials and comparisons use exact integers or rationals.

The program builds each graph from its edge set. It reproduces S(C4) = C8,
the four source edge counts for C4 and P3, and Theorem 41's adjacency,
Laplacian and normalized Laplacian cospectrality for the rook/Shrikhande
pair with K3. The normalized matrix used for the check is I − D⁻¹A,
similar to the source's normalized Laplacian because all checked degrees
are positive.

For the expected positive answer, it checks all 14 labeled regular
graphs of orders 0 through 4: 1,568 output Q matrices, covering every
pair of adjacency spectral classes, all four operations, and both readings
of B. It also checks all four 70-vertex rook/Shrikhande constructions in
both readings, varies the second input separately, and varies both input
pairs in the 128-vertex original-to-original construction.

The survey calls B(G) bipartite but gives edge counts retaining the
original edges, as in the cited R(G). The retained-edge reading is the
primary check; the subdivision reading is checked separately.

These finite checks support fidelity and test the proposed universal
argument. They do not prove the universal assertion. No Lean proof is
included.
