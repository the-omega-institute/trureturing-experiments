# Elliptic generator and bracket checks

`check.py` uses exact SymPy polynomial arithmetic in the quotient
`C[x,y]/(y^2 - (x^3 + a*x + b))`. For the smooth sample curves `(a,b)=(0,1)`
and `(1,1)`, it verifies the source's formula
`[yD,x^nD]=((2n-3)x^(n+2)+(2n-1)a*x^n+2n*b*x^(n-1))D`, the initial seeds and
finite portion of the two-generator recurrence, and the `x^k y D` induction.
It also uses an exact Groebner calculation to rule out a single bracket
`[fD,gD]=D` when both coefficients have x-degree at most two. This is an
independent finite check of the proposed refutation mechanism; it is not a
proof of the uniform claim.

Run with Python 3 and SymPy:

```text
python3 check.py
```

The output is stored in `check.out` and ends with `ALL_CHECKS_OK`.
