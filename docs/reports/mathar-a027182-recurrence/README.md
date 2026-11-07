# OEIS A027182 recurrence check

`check.py` uses exact integer arithmetic for the row self-convolution and the
closed form, and exact SymPy polynomial arithmetic for the variable-coefficient
operator factorization. It reproduces the published terms through `a(5)`,
checks the closed form for `n = 0..25`, and checks the conjectured recurrence
for `n = 5..30`. The OEIS explanatory example prints `120` for `a(2)`, but
`5^2 + 10^2 + 5^2 = 150`, which agrees with the sequence data and definition.

Run from this directory with:

```sh
python3 check.py > check.out
```

The final output line is `ALL_CHECKS_OK` on success.
