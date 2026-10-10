# Wei–Yang Question 6.3: cube polynomial of the associated Mersenne graphs

`check.py N` (SymPy) builds the vertex set M_n of the associated Mersenne graph from the literal definition of arXiv:2407.08237 (circular words in which every run of 1s is immediately followed, circularly, by a strictly longer run of 0s; M_1 = {0}) for n = 1..N. For each n it:
- counts the induced coordinate subcubes Q_k by their top vertex;
- independently counts vertices, edges and induced 4-cycles of the graph itself, without assuming the subcube characterization;
- compares the cube polynomial C_n(x) with the coefficient of z^n in G(x, z) = (z + 2z² + 3xz³ + 5x(1+x)z⁵)/(1 − z − z² − xz³ − x(1+x)z⁵) − 2z²/(1 − z²).

`check-13.log` holds the run for N = 13. Every line ends `OK`, the last line is `ALL_OK`, and the exit code is 0.
