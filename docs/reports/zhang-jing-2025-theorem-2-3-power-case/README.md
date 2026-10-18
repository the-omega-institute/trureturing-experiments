# Zhang–Jing, the open case 1/2 ≤ a ≤ 1 of Theorem 2.3

Check for the open-problem preregistration on the case left open after Theorem 2.3 of F. Zhang and H. Jing, "More inequalities for positive semidefinite matrices", Electronic Journal of Linear Algebra 41 (2025), 595–601: whether [[|A|^{2a}, A*], [A, |A|^{2(1−a)}]] ≥ 0 holds exactly for normal A when 1/2 ≤ a ≤ 1.

Run: `python3 check.py` (Python 3 with numpy).

The program checks, in floating point:
- for a ∈ {0.5, 0.55, 0.6, 0.75, 0.9, 0.99}, n ∈ {2, 3, 4} and ranks n, n − 1, 1: every random nonnormal A gives a block with a negative eigenvalue, and every random normal A gives a positive semidefinite block (1350 pairs);
- at a = 1, the nonnormal A = [[0, 1], [2, 0]] gives a positive semidefinite block, which equals Q*Q with Q = [A, I] (so the block is positive semidefinite for every A);
- at a = 3/4, the singular nonnormal A = [[0, 1], [0, 0]] gives a block with minimum eigenvalue −0.618.

`run.log` holds the output of the run cited by the preregistration.
