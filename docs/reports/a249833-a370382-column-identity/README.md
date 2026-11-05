# A249833 and A370382: exact coefficient check

`check.py` uses Python 3's `fractions.Fraction` and integer arithmetic to
compute the formal solution of `U' = U + U^2 log(U)`, with `U(0) = 1`,
without using the conjectured array. It separately computes the array
`B(0,k)=1` and
`B(n+1,k)=B(n,k+1)+(k+1) sum_{j=0}^k binom(k,j) B(n,j)`.

Sources and attribution:

- Paul D. Hanna, [OEIS A249833](https://oeis.org/A249833/internal), revision
  30 (April 25, 2024): defining EGF and displayed coefficients/expansions.
- Mikhail Kurkov, [OEIS A370382](https://oeis.org/A370382/internal), revision
  15 (November 25, 2024): defining array and downward-antidiagonal entries.
- The column identity is Kurkov's conjecture in both entries, also stated
  in [MathOverflow question 460960](https://mathoverflow.net/q/460960).
- OEIS content is available under its
  [license agreements](https://oeis.org/wiki/The_OEIS_End-User_License_Agreement).

Run from the experiments repository root:

```sh
python3 docs/reports/a249833-a370382-column-identity/check.py
```

The stored `output.txt` reproduces all 23 displayed A249833 terms, its
three displayed related expansions, and all 55 displayed A370382 entries.
It checks `a(0)=1` and `a(n+1)=B(n,0)` for every `0 <= n < 40`, and checks
integrality through degree 40. All checks pass. The final line is
`ALL_CHECKS_OK`.

These are finite consistency and source-fidelity checks. They do not
prove the infinite column identity. The prospective proof uses formal
flow and differential-operator identities; no such proof is implemented
by this program.
