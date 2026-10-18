# Bouthat–Chávez–Sheng, Remark 4.3: the CHS norm of the path

Checks for the open-problem preregistrations on Remark 4.3 of L. Bouthat, Á. Chávez and S. Sheng, "A general framework for inequalities on simple graphs", arXiv:2606.17356v1.

Run: `python3 check.py 30 60` (Python 3 with mpmath).

The program computes ‖P_n‖_{X,d}^d = h_d(λ(P_n)) = [t^d] 1/det(I − tA_n) exactly, using det(I − tA_n) = det(I − tA_{n−1}) − t² det(I − tA_{n−2}). For 1 ≤ n ≤ 30 and even 2 ≤ d ≤ 60 it compares this value with four expressions:
- the trigonometric expression, at 80 digits;
- the bilateral binomial sum, exactly;
- the stated truncation of the bilateral sum;
- the closing closed form under n + 1 ≥ d/2.

`run.log` holds the output of the run cited by the preregistrations.
