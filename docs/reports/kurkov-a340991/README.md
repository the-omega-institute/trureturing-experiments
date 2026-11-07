# OEIS A340991 exact preregistration check

`check.py` uses Python integers only. It independently computes the finite Cauchy-product coefficients and the shifted Kurkov array, then checks their polynomial identity for the prime sequence and for a signed arbitrary sequence. It also reproduces OEIS row 4, `[0, 7, 29, 36, 16]`. The stored output is `check.out`.
