# Bardakov–Kozlovskaya Conjecture 3.3(1)

`check.py` is an exact SymPy check for arXiv:2305.07261v1, Conjecture 3.3(1). It constructs the source polynomial `p_n` from the root product by an integer resultant, rewrites it in the elementary symmetric variables, reproduces the source tables for `n = 4, 6, 7`, and checks `p_(p^m) ≡ e1^(p^m) (mod p)` for `n = 2, 3, 4, 5, 7, 8, 9`. The `n = 6` control confirms that the scope is prime powers. All arithmetic is exact integer/SymPy arithmetic; no floating point is used.

Run:

```sh
python3 check.py > check.out
```

The final output line is `ALL_CHECKS_OK`.
