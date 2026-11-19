# Balabdaoui–Wellner Gaussian factor

Run `python3 check.py` with Python 3.12 or later and SymPy (verified with SymPy 1.14.0). Expected output is in [check.out](check.out).

The fidelity check reproduces the Gaussian product in equation (4.1) through its logarithmic derivative and its normalization at zero. It takes the derivative of the integral in (4.2) by the fundamental theorem of calculus. The settlement check generates the first four Taylor coefficients by a recurrence and by direct series expansion. It verifies the strictly positive polynomial that contradicts the residue identity forced by the Schoenberg product. Analytic continuation, the Schoenberg implication, and comparison of residues enter the mathematical argument.

Source is Balabdaoui and Wellner, *Chernoff’s density is log-concave*, §4(b), arXiv:1203.0828v2. [Paper](https://arxiv.org/abs/1203.0828v2)
