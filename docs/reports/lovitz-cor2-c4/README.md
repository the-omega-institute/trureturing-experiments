# Lovitz `Cor_2(C^4)` exact check

`check.py` uses only SymPy exact arithmetic. It reproduces the source convention that a correlation matrix is positive semidefinite with unit diagonal on the 2 by 2 identity, then checks a nontrivial rank-three four-vector family: four explicitly normalized product vectors in `C^2 ⊗ C^2` give a rank-three Gram matrix `P`, and the two local Gram matrices `Q` and `R` have rank at most two with `P = Q ⊙ R` entrywise exactly.

Run:

```sh
python3 check.py
```

The final line is `ALL_CHECKS_OK`.
