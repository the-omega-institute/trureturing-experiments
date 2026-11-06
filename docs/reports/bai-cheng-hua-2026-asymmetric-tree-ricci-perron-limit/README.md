# Asymmetric tree Ricci subdivision checks

Reference: Shuliang Bai, Haoxuan Cheng and Bobo Hua, *Edge Subdivision and the Perron Eigenvalue of Tree Ricci Matrices*, arXiv:2605.30949v1, Remark 5(2) and Section 5.2. Source version submitted 29 May 2026.

Run `python3 check.py` with Python 3 and SymPy. The script uses integers, exact rational arithmetic and exact polynomial algebra throughout. It prints `ALL_CHECKS_OK` on success. The output is in `output.txt`; all 600 rational spectral intervals are in `intervals.csv`.

The direct Ricci-matrix construction reproduces the positive zero eigenvector and interface data of the source's S_(2,5) example. Its characteristic polynomial agrees with an independently constructed leaf-orbit quotient for three families at lengths 1 through 4. A diagonal similarity transforms this quotient to a real symmetric Jacobi matrix with the stated off-diagonal products. Sylvester's criterion and rational bisection then enclose its largest eigenvalue, which is the original matrix's Perron eigenvalue. The remaining leaf-difference eigenvalues are negative.

For each of S_(2,5), S_(4,6), and S_(1,4), lengths 1 through 200 are checked. The first and third families increase, while S_(4,6) decreases. Every length-200 eigenvalue is within 1/10^10 of its nominated algebraic limit. The S_(2,5) limit satisfies 216 L^3 + 288 L^2 + 111 L - 2 = 0 and lies strictly between 1/60 and 1/50. Its branch equation is identical to the S_(5,5) equation. The flat S_(3,3) sanity check covers lengths 1 through 200.

These are finite checks supporting a proposed general argument. They do not prove convergence for every tree or substitute for a Lean settlement. The source's universal conjecture assumes both positive matching roots exist uniquely; the S_(2,5) example has only one positive matching root.
