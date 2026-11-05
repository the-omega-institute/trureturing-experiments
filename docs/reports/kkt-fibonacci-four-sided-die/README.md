# Four-sided die and Fibonacci targets

`check.py` uses only exact integers, rational fractions, and SymPy’s exact
quadratic expression for the golden ratio. It reproduces the source’s
four-sided heuristic `(3/5)φ < 1`, verifies the renewal masses and the
`5/13` minorant used by the proposed proof, checks the Fibonacci potential
telescope, and independently enumerates finite stopped prefixes to verify the
geometric survival contraction for the first tested thresholds.

Run:

```sh
python3 check.py
```

The final line is `ALL_CHECKS_OK`.
