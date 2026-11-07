# Yoshida Open Problem 5: rapid decrease under `S_r`

`check.py` uses only exact integers, `fractions.Fraction`, and SymPy exact
algebra. It reproduces the 4-rapid-decrease identity in Yoshida's Example 30
for `a_k = 2^(-k^2)` and checks the refuting family

`a_k = L (s/L)^k A^(-k(k-1)/2)`, `s=q^2-1`, `A=(q+1)^2/2`.

For `q=7`, `L=49`, and `r=1`, the output's exact minimum adjacent ratio is
`A^2(1-A^(-1))`, while the source threshold is strictly larger than `b_1/2`;
the polynomial certificate in `u=q-7` is also checked symbolically. Run:

```sh
python3 check.py
```

The final line is `ALL_CHECKS_OK`.
