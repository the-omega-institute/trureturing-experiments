# OEIS A393700 odd-value multiplicity check

`check.py` uses exact integer arithmetic. It computes the strict inside/outside unit-square count from the four corner squared norms, reproduces the first 21 displayed terms of OEIS A393562, verifies the fourfold reflection reduction to the nonnegative quadrant on small radii, and checks the proposed two-witness construction for every `s` from 2 through 250. The construction yields the target odd value `2*s - 1` at two distinct radii. No floating point arithmetic is used.

Run with:

```text
python3 check.py
```

The stored output is `check.out`; the final line is `ALL_CHECKS_OK`.
