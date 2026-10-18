# Cullinan–Young (arXiv:2603.05408), Remark 8: the leading coefficient of I_N

Source: "Numerical experimentation suggests that the leading coefficient of I_N(x) (and, therefore, F_N(x)) is given by (−1)^{(N−2)/2} C((N−2)/2)/(N−1)!, where C(n) is the nth Catalan number." Here I_N(x) = Σ_{i=0}^{N} sgn(i − N/2) ℓ_i(x) for even N, and ℓ_i are the Lagrange basis polynomials on the nodes j − N/2.

Route (N = 2m):
- **Degree.** The x^{2m} coefficient cancels in symmetric pairs i ↔ 2m − i.
- **The next coefficient.** Because the nodes sum to 0, the x^{2m−1} coefficient of each Lagrange numerator is its omitted node a_i. Pairing reduces L = [x^{2m−1}] I_{2m} to 2S/(2m)!, where S = Σ_{k<m} (m − k)(−1)^k C(2m, k).
- **The sum.** Two telescoping alternating-binomial identities give S = (−1)^{m−1} C(2m − 2, m − 1). Hence L = (−1)^{m−1} C_{m−1}/(2m − 1)!.

## Program

`check.py` (standard library only). Run `python3 check.py`. It expands I_N exactly for N = 2, 4, …, 40, checks that deg I_N = N − 1 and that the x^{N−1} coefficient equals the formula, and prints `RESULT PASS`.

Reading: `run.log` (exit 0, `RESULT PASS`). The first values are 1, −1/6, 1/60, −1/1008.
