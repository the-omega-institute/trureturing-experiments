# Complex symmetric big-centroid interpolation check

`check.py` verifies the source construction for `n = 2, 3` and `d = 3, 4` using exact rational polynomial arithmetic in SymPy.  It finds a rational `c` for which

\[
P(z)=Q(z)(z^{d-2}+c)+Q'(z)
\]

is squarefree and coprime to `Q`, checks `P(z_k)=Q'(z_k)`, and verifies the Lemma 3 coefficient identity

\[
q_I(\partial)\sum_k x_k^d=d!\,\ell_I
\]

through the exact factorization `P(z_k)=H_I(z_k)\prod_{i\notin I}(z_k-t_i)`.  The same exact conditions establish the off-diagonal interpolation zeros and nonzero diagonal factors.  No floating-point value is used as evidence.

Run:

```sh
python3 check.py
```

The captured output is in `check.out`; a successful run ends with `ALL_CHECKS_OK`.
