# Duplication-corona normalized spectra

Exact finite diagnostics for Sason, Krupnik, Hamud and Berman, *On Spectral Graph Determination*, arXiv:2412.20775v6, Question 5.19 (journal Question 2, DOI 10.3390/math13040549).

Run `python3 check.py` with Python 3 and SymPy 1.14.0. The program uses only integers, exact symbolic algebra and `fractions.Fraction`; `output.txt` contains its output.

The source convention is D^(-1/2)(D-A)D^(-1/2), with inverse square root zero at an isolated vertex. The program computes the characteristic polynomial through the similar rational matrix D^(-1)(D-A) on positive-degree vertices, splitting off zero rows for isolated vertices. It clears denominators before computing each integer characteristic polynomial and rescales its coefficients exactly. A direct symbolic square-root calculation checks the convention on a graph with isolates.

Fidelity: reproduce the A, L and Q spectral equalities of Theorem 5.18 / journal Theorem 39 with G=G'=K2 and the nonisomorphic rook/Shrikhande pair as attached factors. The same three constructions have equal normalized characteristic polynomials. Further tests vary the base factors (orders 64, 64, 128) and both factors at once (orders 120, 120, 220), and exhaust all regular labeled graphs with 0–4 vertices: 14 graphs, 10 adjacency-spectrum groups, covering all 676 ordered factor quadruples by 588 common-reference comparisons. The small cases check the stronger identity without the source's attached-factor nonisomorphism restriction.

The rook graph has vertices in Z4×Z4 adjacent in a common row or column; Shrikhande uses differences ±(1,0), ±(0,1), ±(1,1). Their common adjacency polynomial and distinct K4 counts are checked. The 10-vertex pair uses the mathematical neighborhood table in Igal Sason, arXiv:2609.10884, also illustrated in the source survey's Figure 1; its cospectrality and distinct edge-triangle invariants are checked independently. No third-party program code is copied.

These are finite diagnostics, not a proof of the universally quantified question and not Lean verification. The general argument must establish compatible vertex/unsigned-incidence transports and the three class-constant degree computations, including degree zero.
