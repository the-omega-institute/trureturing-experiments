# Sun's Tsallis heat-flow lower endpoint

The target is the open lower-endpoint question in Section VI and Remark 20 of Lukang Sun, *Concavity of Tsallis Entropy and Tsallis Entropy Power along Heat Flow*, arXiv:2601.12944v3. The proposed Euclidean threshold is

\[
q_d=\frac{2(d-1)^2}{2d^2+1},\qquad d\ge2.
\]

Run `python3 check.py` with Python 3 and SymPy 1.14.0. The stored output is `output.txt`. Every comparison uses exact integers, fractions or symbolic algebra; no floating point or numerical quadrature supplies evidence.

The fidelity check reproduces Sun's Remark 20 formula for the entropy's second derivative on the Gaussian heat kernel, with generator Δ. The conditional algebra checks the trace-free square identity and the threshold factorization. For each integer dimension 2 through 12, the program verifies the radial Laplacian and equality condition for u(x)=(1+|x|²)^(-(d-1)), evaluates D and A through exact beta integrals, and checks A/D=-(d-1)²/[d(d+2)]. At three rational fractions of each threshold it checks normalizability, finite initial entropy and the strict initial sign. The positive angular factor omitted from D and A is |S^(d-1)|/2.

The square identity and the associated exponent are classical ingredients: Dolbeault–Esteban–Loss, arXiv:1302.0814, Lemma 10, and Cai–Wang, arXiv:2512.15525v2, Theorem 4, must be credited. This computation does not establish literature novelty, the analytic integration by parts for general densities, differentiation under heat convolution, or persistence of the initial strict sign at positive times. It is a preregistration check, with no Lean proof.
