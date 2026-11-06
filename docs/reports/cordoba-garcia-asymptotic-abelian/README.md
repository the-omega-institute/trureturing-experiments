# Cordoba-Caycedo--Garcia, Conjecture 4.30

`check.py` uses exact integer and rational arithmetic to verify the primary-decomposition formulas for finite abelian groups, the equality `|G[n]| = |G/nG|`, and the bounded-width obstruction used in the preregistered proof route.  The family `G_r = direct_sum_{j=1}^r Z/2^j Z` has 2-width one but unbounded 2-rank; for `N = 2` its squared torsion-to-scale ratio is computed exactly and tends below the required positive asymptotic scale.  It also checks the fixed-kernel bound for all cyclic groups `Z/m` with `m < 80` and `n < 12`.

Run:

```sh
python3 check.py
```

The final output line is `ALL_CHECKS_OK`.
