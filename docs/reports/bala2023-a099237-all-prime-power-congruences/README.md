# Bala's A099237 prime-power conjectures

Source: [OEIS A099237](https://oeis.org/A099237), published revision #21,
February 3, 2026. Paul Barry defined the sequence in 2004; Peter Bala's
three conjectures were contributed January 19, 2023.

The sequence is `a(n) = sum(C(n*(n-j), j), j=0..n)`, using ordinary
nonnegative-integer binomial coefficients. The proposed congruences are
`a(2^k) = 1 mod 2^k`, `a(3^k) = 1 mod 3^(k+1)`, and
`a(p^k) = 1 mod p^(k+1)` for every prime `p >= 5`.

Run with Python 3.12 or newer, from this directory:

```sh
python3 check.py > check.out
```

The standard-library-only checker uses exact integers throughout. It reproduces
all 23 terms printed in the source, checks the equivalent reindexed formula
through index 200, and verifies all 60 positive-exponent prime-power cases with
`p^k <= 200` (5025 interior binomial terms). It separately checks the endpoint
`k=0`, the Frobenius coefficients, the alternating binomial row, the absorption
identity, valuation bookkeeping, and cancellation of the unit factor. It also
checks the stronger dyadic helper and the failure of a uniform extra odd-prime
power at `p=7, k=1`.

`check.out` is the deterministic output. Its final line is `ALL_CHECKS_OK`.
Finite checks support the preregistered general argument; they do not prove its
unbounded conclusion. No Lean proof is included.

The checker is independently authored and licensed Apache-2.0. Source terms and
attribution are from the OEIS under its [end-user licence](https://oeis.org/LICENSE).
