# Constant Wronskians with prescribed finite poles

Run `python3 check.py` with Python 3.12+ and SymPy 1.14.0. The stored output is in `check.out`.

The program computes the Section 2 polynomial example from Hermoso and Alcázar and verifies W = −8. It then differentiates explicit rational-function columns, clears their row denominators and computes determinants over the exact rational polynomial ring. For each family member it checks the reduced denominator, independence by distinct numerator degrees over a common denominator and the predicted nonzero factorial product. The tests include dimensions 2 through 5, both integer and nonintegral rational pole sets and a case that requires a nonzero shift of the antiderivative. These are finite checks of the construction in the preregistration. [arXiv:2410.18867v1](https://arxiv.org/abs/2410.18867v1)
