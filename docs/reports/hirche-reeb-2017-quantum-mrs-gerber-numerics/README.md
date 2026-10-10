# Hirche–Reeb Conjecture VII.1 (quantum Mrs. Gerber's lemma): qubit numerical search

`opt.py seed d restarts` minimizes H(X₁⊕X₂|B₁B₂) minus the conjectured two-branch lower bound of arXiv:1706.09752v2, eq. (72). It works over independent classical–quantum inputs with arbitrary priors and arbitrary d×d output density matrices, using natural logarithms (NumPy, SciPy 1.18.1; Nelder–Mead then BFGS).

`opt.py 1 2 30` (qubits, 30 restarts) gives a minimum gap of −8.4e-12, which is numerical round-off at the equality cases. No counterexample was found in this scope.
