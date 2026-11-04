# Duhamel–Cazenave logarithmic adjacency bound

Source: Elie Duhamel and Tristan Cazenave, *Spectral Graph Theory Conjecture
Generation*, AAAI 2026 Workshop AI for Research,
[author-hosted PDF](https://www.lamsade.dauphine.fr/~cazenave/papers/SpectralGraphTheoryConjectureGeneration.pdf),
Table 3, row with runtime 45 and score −0.15645678373681005:
the second largest adjacency eigenvalue is conjectured to be at most twice the
natural logarithm of the largest combinatorial Laplacian eigenvalue.
The paper concerns connected undirected graphs and defines the Laplacian as D−A.

Run `python3 check.py` with Python 3 and SymPy. The stored output is `output.txt`.
The checked environment uses Python 3.14.8 and SymPy 1.14.0.

The fidelity check reproduces source equation (1),
λmax(A) ≤ λmax(L)−1, with equality on complete graphs of orders 2 through 10.
It computes their exact characteristic polynomials. It does not reproduce the
source's randomized search, runtime or floating point score.

For Gm = Km □ K2, the program independently builds adjacency and degree-minus-
adjacency matrices and checks a full integer eigenbasis, both characteristic
polynomials, connectivity, degrees and edge counts for 2 ≤ m ≤ 14. The adjacency
spectrum is m, m−2, 0^(m−1), (−2)^(m−1), and the Laplacian spectrum is
0, 2, m^(m−1), (m+2)^(m−1). It also checks the symbolic subspace-action formulas.
The complete-graph spectra and Cartesian-product spectral rule are classical;
see Brouwer–Haemers, *Spectra of Graphs*, §§1.4.1 and 1.4.6.

The exact exponential certificate 443/48 > 9 gives 2 ln9 < 5. Together with
ln(1+x) ≤ x, the symbolic slack identity
(m−2)−(5+2(m−7)/9) = 7(m−7)/9 proves the proposed comparison for every integer
m ≥ 7. At m=7, the full basis determinant is −6272, the graph has 14 vertices
and 49 edges, and 5 > 2 ln9 contradicts the selected conjecture.

An additional tree check verifies quotient matrices and all leaf-difference
eigenvector equations for two joined centers with 256 leaves each. It records
exact bounds giving a second adjacency eigenvalue greater than 31/2 and twice
the logarithm of the largest Laplacian eigenvalue less than 14.

All evidence uses integers, rational numbers and symbolic polynomial algebra.
The analytic implications use the positive exponential series and ln(1+x) ≤ x;
the program does not implement a theorem prover for those analytic facts.
These checks support a preregistration outline and are not a Lean kernel proof.
