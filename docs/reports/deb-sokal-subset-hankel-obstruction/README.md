# Deb–Sokal subset-Hankel obstruction

`check.py` uses exact integer arithmetic. It implements the recurrence for the higher-order Stirling subset numbers
\[
S^{(r)}_{n,k}=\binom{n+(r-1)k-1}{r-1}S^{(r)}_{n-1,k-1}+kS^{(r)}_{n-1,k},
\]
with `S^(r)_(0,0)=1`, forms the row polynomials `s_(r,n)(x)`, and evaluates the 3×3 Hankel determinant for rows `m,m+1,m+2` and columns `0,1,2`. It reproduces the source’s Appendix C witness (`r=3`, `m=1`) and the tabulated `r=4,n=2` row, then checks the all-`r`, all-`m` coefficient obstruction and the finite-difference identities for `r=3..15`, `m=1..15`. No floating point is used.

Run:

```sh
python3 check.py > check.out
```

The final line is `ALL_CHECKS_OK`.
