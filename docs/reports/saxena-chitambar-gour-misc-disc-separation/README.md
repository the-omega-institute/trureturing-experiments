# Exact coherence-distillation checks

Run `python3 check.py` with Python 3 and SymPy. The saved output is `check.out`.
The program uses exact integer, rational, and algebraic arithmetic throughout.

It reproduces the equal optimal fidelities for the partially depolarizing and
partially dephasing channels in Saxena, Chitambar and Gour,
*Dynamical Resource Theory of Quantum Coherence*, arXiv:1910.00708v4, Section V.E,
DOI: 10.1103/PhysRevResearch.2.023298, for 128 parameter triples.
The eigenvalue constraints determine an exact feasible interval, and its
endpoints independently certify the optimum of the reduced linear objective.

For the Fourier measure-and-prepare channel
Φ_c(A) = Σ_j Tr(P_j A)(I + c_j X)/2, the program checks 48 response profiles in
dimensions 2, 3, 4 and 6. It verifies Fourier projectors, the input-first Choi
matrix and its trace-preserving normalization, exact spectra, feasible MISC and
DISC matrices, cyclic invariance, and the optimum of the reduced off-diagonal
box. A symbolic positivity identity certifies the off-diagonal bound; symbolic
identities give the dimension-uniform separation margins for c = (q,0,…,0),
1/2 < q < 1. The strict threshold is the defining equation's > 1 − ε.

These are certificates and finite checks supporting a preregistration outline.
They are not a Lean proof of the operational SDP equivalence, the general
cyclic reduction, or the source's open problem. The planned proof is a general
argument in arbitrary input dimension, rather than positive finite enumeration.

Source TeX: https://arxiv.org/e-print/1910.00708v4
Journal: https://doi.org/10.1103/PhysRevResearch.2.023298
No third-party implementation is copied.
