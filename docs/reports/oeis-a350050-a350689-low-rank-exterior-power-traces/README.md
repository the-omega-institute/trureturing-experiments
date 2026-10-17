# OEIS A350689 and A350050: Spezia's exterior-power trace conjectures for three low-rank matrices

S. Spezia conjectures that the traces of exterior powers of three structured n × n matrices are given by A350689 and A350050 and vanish beyond a fixed order. The entries define these traces as, up to sign, the coefficients of the characteristic polynomial. The three matrices are:
- **A350689.** M[i,j] = ⌊(j − i + 1)/2⌋. The conjecture: |trace Λ³M| = A350689(n) and trace Λ^k M = 0 for k > 3.
- **A350050, matrix of A349107.** A[i,j] = n − |(n+1)/2 − j| − |(n+1)/2 − i|. The conjecture: |trace Λ²A| = A350050(n) and trace Λ^k A = 0 for k > 2.
- **A350050, matrix of A349108.** B[i,j] = (n mod 2) + |(n+1)/2 − j| + |(n+1)/2 − i|. The same two conjectures.

The proof is that the matrices have rank ≤ 3, respectively ≤ 2:
- M[i,j] = (j − i)/2 + 1/4 − (−1)^(i+j)/4;
- A and B are sums of a constant matrix and two rank-one matrices built from g_i = |(n+1)/2 − i|.

With M = UVᵀ, charpoly(UVᵀ) = x^(n−r)·charpoly(VᵀU). This leaves a 3 × 3 (or 2 × 2) determinant in closed form.

`check.py` (Python 3 with sympy) verifies for n ≤ N:
- the factorizations;
- the exact characteristic polynomials;
- the closed forms of det(VᵀU);
- the published data.

Run `python3 check.py 22`. Exit 0 and `RESULT PASS` mean every check holds.
