# Sampling overhead and Schur convexity

Exact checks for the unnumbered conjecture following Eq. (53), immediately before Proposition 4, in Xiong, Chandra, Ng and Hanzo, *Sampling Overhead Analysis of Quantum Error Mitigation: Uncoded vs. Coded Systems*, arXiv:2012.08378v1; DOI 10.1109/ACCESS.2020.3045016.

Run `python3 docs/reports/sampling-overhead-schur-convexity/check.py` from the repository root. Python 3 and SymPy are required. The output is stored in `output.txt`. Inputs, algebra, inversions and comparisons use integers, rational numbers and exact symbolic expressions.

The fidelity check reproduces the source's single-qubit depolarizing example: Pauli transfer eigenvalues 1 and 24/25, inverse coefficient norm 17/16, and sampling overhead factor 33/256. It checks all 256 two-qubit Pauli multiplication rules modulo phase against XOR labels.

The family check verifies the symbolic inverse and Fourier sign-test identities and checks 30 rational parameter choices in dimensions 16 and 64. The accompanying general argument uses 0 < a < 1/6 and 0 < t < 1, with uniform-mixture parameter d = a²t(1−6a). The rational grid is supplementary evidence, not a verification of all real parameters.

The interior witness has probability numerators (5549,797,797,797,5,…,5) and (5549,797,797,401,401,5,…,5), both with denominator 8000. Averaging labels 3 and 4 decreases the fourth descending prefix sum by 99/2000, preserves generalized gate error probability 2451/8000, and increases the sampling overhead factor by 22957327/50808384. Fourier inversion and independent rational Gaussian elimination agree, and all convolution equations are checked.

These checks support a proposed refutation of the unrestricted conjecture. They do not constitute a Lean proof or refute the local tensor-product statement of Proposition 4.
