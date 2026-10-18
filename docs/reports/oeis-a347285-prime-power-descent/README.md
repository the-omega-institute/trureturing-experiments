# OEIS A347285: descending prime-power exponents

Check for the open-problem preregistration on De Vlieger's conjecture in OEIS A347285: for n > 0, n − S = 1, where S is the sum of the absolute first differences of row n.

Run: `python3 check.py 400` (Python 3, standard library only).

Row n starts with n. Each next term is e_k = ⌊log_{p_k}(p_{k−1}^{e_{k−1}})⌋, computed exactly as the largest e with p_k^e ≤ p_{k−1}^{e_{k−1}}, and the row continues while e_k > 0.

For rows 1, …, 400 the program checks:
- the published data;
- that each row is strictly decreasing and ends in 1;
- that n − S = 1.

`run.log` holds the output of the run cited by the preregistration.
