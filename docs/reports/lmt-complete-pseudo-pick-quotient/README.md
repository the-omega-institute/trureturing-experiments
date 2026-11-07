# Complete pseudo-Pick pairs and quotient positivity

The proposed counterfamily concerns Luo, McCullough and Tsikalas, *Beurling Criteria for Reproducing Kernels*, arXiv:2607.25088v2, Question 9.3. The parameter is a real number satisfying 0 < epsilon < (sqrt(19) - 3)/10. The core consists of epsilon, 1, 2, with countably many additional labels indexed by positive integers.

Run `python3 check.py` with Python 3 and SymPy 1.14.0. The stored output is `output.txt`. The script uses integers, rational arithmetic and exact symbolic algebra; it uses no floating point arithmetic as evidence.

The fidelity check reproduces the exact rational decomposition and first three Taylor coefficients in source Example 8.3 (TeX label `eg:recursivel`). The remaining checks verify the symbolic negative quotient form, the core rigidity relations, unrestricted complex symbolic block identities at sizes 1 and 2, anchor-compression coefficient identities, and 892 exact finite quotient submatrices. Extra labels include nonconsecutive integers.

These computations support the proposed counterexample and detect errors in its algebra. They do not establish interpolation for all finite node sets and all matrix sizes. That conclusion requires the general null-form, Schur-complement and compression argument in the preregistration, together with the classical complete Pick theorem. No Lean proof is included.

Source: https://arxiv.org/abs/2607.25088v2 . Complete Pick ingredient: Agler and McCarthy, *Complete Nevanlinna-Pick Kernels*, J. Funct. Anal. 175 (2000), 111–124, Theorems 3.1 and 4.2, https://doi.org/10.1006/jfan.2000.3599 . The program is independently written and contains no copied third-party code.
