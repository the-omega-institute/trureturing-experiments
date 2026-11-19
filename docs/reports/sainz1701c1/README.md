# Sainz–Wolfe Conjecture 1

Run `python3 check.py` with Python 3.12 or later and SymPy. The saved output in `check.out` was produced with SymPy 1.14.0.

The fidelity check transcribes the two incidence tables in Appendix E of arXiv:1701.05171v2 and independently generates the 12-context bipartite Foulis–Randall product. It compares augmented row spaces and tests each signed row dependence exactly. The family check tests the added graph edges, incidence identity and symbolic projection certificate for every block-size tuple in `{1,2,3}⁴`, including both zero-mass cases. The final check constructs the 20-event quantum model and checks its probabilities and CHSH violation with exact algebra.

The equality of the entire normalized model sets for arbitrary cores uses the tensor construction and certificate restriction in the preregistration argument. The program checks specific witnesses and identities.
