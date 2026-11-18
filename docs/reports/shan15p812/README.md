# Shan’s inverse questions

Problem 8.1.2 of Shen Shan’s 2015 Western University dissertation asks whether the compositional inverse of `f(x)/x` is completely monotone for completely monotone f, and whether the compositional inverse of `xg(x)` is Bernstein for Bernstein g. [Thesis](https://hdl.handle.net/20.500.14721/36017)

Run `python3 check.py` with Python 3.12 or newer and SymPy 1.14.0. The script reproduces the derivative recurrence in thesis equation (2.2) for a two-atom measure. It uses exact symbolic differentiation to obtain the inverse formulas, checks the positive third inverse derivative for the two-atom family at every real T ≥ 30, and checks the positive second inverse derivative for the atomic Lévy family throughout the proposed drift interval. The output is in `check.out`. The checks cover the algebra and sign calculations. The global inverse arguments are given in the preregistration draft.
