# Two-conformal rigidity without one-conformal rigidity

Source: H. Assumpção, G. Coutinho and C. Godsil, *Total Conformal Rigidity in Graphs*, arXiv:2605.08508v2, Definition 1.2 and the unnumbered question in Section 9, printed p.20.

Run `python3 check.py` with Python 3 and SymPy 1.14.0. All evidence uses exact integers, rationals or symbolic polynomials. The stored output ends in `ALL_CHECKS_OK`.

The fidelity check reproduces equations (1) and (2): the weighted Laplacian, its adjoint edge quadratic form, and the symbolic trace pairing for arbitrary weights and an arbitrary symmetric test matrix. It also checks zero row sums and twice-the-total-weight trace.

For the proposed graph K₅ minus edge {1,2}, the program checks the uniform spectrum (0,3,5,5,5), the rational matrix X satisfying 0 ≤ X ≤ I and tr X = 2, and the identity tr(L(w)X) = (10/9)∑w. The identity is symbolic in all nine weights, so it covers every point of the real weight simplex, including its boundary. It also checks Y = I − J/5 − X, a feasible lower-sum certificate with tr(L(w)Y) = (8/9)∑w.

The general variational argument required to interpret these certificates is: for any symmetric matrix A and positive contraction X of trace k, tr(AX) is at most the sum of its k largest eigenvalues. In an orthonormal eigenbasis, the diagonal entries of X lie in [0,1] and sum to k; allocating this mass to the largest eigenvalues gives the bound. Thus S₂(w) ≥ 10 for every admissible w, and s₂(w) = 18 − S₂(w) ≤ 8. Uniform weights attain both. This universal argument is not a finite enumeration.

The positive weights 6/5 on the six edges meeting {1,2}, and 3/5 on the remaining three edges, have total 9 and spectrum (0,18/5,21/5,21/5,6). Their second eigenvalue is strictly greater than the uniform value 3. The program additionally checks all nine vertices of the weight simplex, where the weighted graph is disconnected.

These are exact certificate checks supporting a preregistration outline. No Lean proof of the settlement is included. The infinite-family generalization is outside this report's target.
