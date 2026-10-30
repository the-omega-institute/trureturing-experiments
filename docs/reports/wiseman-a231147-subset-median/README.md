# OEIS A231147 subset-median check

`check.py` uses only exact integer arithmetic. It expands
\(N_n(x)=\sum_{j=0}^{n-1}x^{n-1-j}(1+x^2)^j\), enumerates every nonempty subset of
\(\{1,\ldots,n\}\), and compares the coefficient row with the doubled-median
index for every `n` from 1 through 16. It also checks the complement-and-shift
encoding and inverse for every subset and every pair `(j,Q)` in that range.

Run:

```text
python3 check.py
```

The output ends with `ALL_CHECKS_OK` and reports 131054 checked nonempty
subsets. The `n=5` row is `[1, 1, 5, 4, 9, 4, 5, 1, 1]`, matching the OEIS
example.
