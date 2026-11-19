# Gy’s general coefficient log-concavity question

René Gy asks whether the two triple-product differences defined in Conjecture 5.1 have log-concave coefficient sequences for arbitrary strongly q-log-concave inputs. The question appears in §5, printed page 16. [Integers 26 (2026), A108](https://math.colgate.edu/~integers/aa108/aa108.pdf)

Run `python3 check.py` with Python 3.12+ and SymPy 1.14.0. The check enumerates the source’s symmetric monomials, reproduces the natural q-Stirling specialization against independent recurrences, and checks the adjacent-triple counterexample using exact identities and direct rational instances. The default run output is in `check.out`.
