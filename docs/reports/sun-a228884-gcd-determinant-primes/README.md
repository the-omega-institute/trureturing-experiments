# Sun's OEIS A228884: composite-index determinant prime divisors

Source: [OEIS A228884](https://oeis.org/A228884/internal), revision #17.
The matrix is `C_n(i,j) = gcd(i-j,n)` for `1 <= i,j <= n`, and `D_n = det(C_n)`.
Clause (ii) conjectures that every prime divisor of `D_n` is smaller than a
composite index `n`. This report concerns only clause (ii).

Run with Python 3, using only its standard library:

```sh
python3 docs/reports/sun-a228884-gcd-determinant-primes/check.py
```

`output.txt` contains the exact output. The checker uses integer Bareiss
elimination on the actual matrices for `0 <= n <= 60`, with independent rational
Gaussian elimination through `n = 12`. It reproduces all 21 OEIS DATA terms
`n = 0..20` and the source's prime-index identity
`D_p = (p-1)^(p-1) (2p-1)` through 60.

For every `1 <= n <= 60`, the actual determinant agrees with both the integer
Fourier eigenvalue product and the prime-power determinant product. All composite
indices in that range satisfy the conjectured prime bound. The program also
checks every CRT entry identity and determinant identity for ordered coprime
pairs `uv <= 60`, and the local zero-frequency bound for primes `p <= 97`,
`2 <= a <= 10`. Every Bareiss division and every rational denominator is checked.
The final output line is `ALL_CHECKS_OK`.

These are finite checks, not a proof of the universal assertion. The Fourier
formulas are existing ingredients: van der Kamp, arXiv:1201.3139v1, equation (6)
and the prime-power formula (16); Holleboom, arXiv:1410.2054v1, Theorem 1. No
floating-point computation or random sampling is used.
