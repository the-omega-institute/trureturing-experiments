# Zhang–Jing, Question 2.14: the functional block-normality criterion

Check for the open-problem preregistration refuting Question 2.14 of F. Zhang and H. Jing, "More inequalities for positive semidefinite matrices", Electronic Journal of Linear Algebra 41 (2025), 595–601.

Run: `python3 check.py` (Python 3 with numpy; the s = 2 case uses exact rationals).

With g(x) = 1 + (x − 1)/(1 + (x − 1)²), f(x) = x²/g(x) and A_s = [[0, 1 + 1/s], [1 + s, 0]], the program checks:
- the admissibility of f and g on a grid;
- the exact case s = 2: |A| = diag(3, 3/2), g(3) = g(3/2) = 7/5, c·M = QᵀQ with Q = [A, cI], and A nonnormal;
- the family for several s: M ≥ 0 and A nonnormal;
- a contrast with an injective admissible pair (f = x^{3/2}, g = x^{1/2}), where M is not positive semidefinite.

`run.log` holds the output of the run cited by the preregistration.
