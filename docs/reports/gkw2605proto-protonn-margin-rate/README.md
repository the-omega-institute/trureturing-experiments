# Proto-NN margin-rate preregistration check

`check.py` uses only exact integers and `fractions.Fraction`. It reproduces the witness distribution's Bayes risk `L*=1/4`, the three-prototype occupancy expectation, the one-label-cell regret bound, and the lower bound `(k/n)^3/16384` for `k≥4`, `n≥4k`. For `gamma=7` and `n=t^5`, `k=t^4`, the conjectured terms are `4/t^4` and `1/t^4`, while the lower bound is `1/(16384 t^3)`, so their ratio grows as `t/81920`.

Run:

```sh
python3 check.py
```

The expected final line is `ALL_CHECKS_OK`; the captured output is in `check.out`.
