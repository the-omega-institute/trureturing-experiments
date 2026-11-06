# OEIS A396563: Stirling minor and nested sum

Source: Mikhail Kurkov, [OEIS A396563](https://oeis.org/A396563/internal),
approved revision 15, stamp 1 June 2026 13:43:04. The conjecture is the sole
FORMULA line. The numeric fidelity data are the six displayed EXAMPLE rows;
OEIS data retain the [OEIS end-user licence](https://oeis.org/LICENSE) and
attribution. The verification program is independently written.

Run with Python 3; no third-party packages are needed:

```sh
python3 docs/reports/oeis-a396563-stirling-minor/check.py
```

The program uses integers and `fractions.Fraction` throughout. It reconstructs
unsigned Stirling numbers from their recurrence and determinants by rational
Gaussian elimination. The conjectured side is evaluated literally by recursive
nested-chain enumeration, falling products, and factorial division.

Results in `check.out`:

- All 42 displayed source values are reproduced, including B(2,2,2) = 725.
- The literal conjecture holds on all 245 triples with 0 <= n,k <= 6 and
  1 <= m <= 5, including empty determinants and n = 0.
- Independent subset enumeration, elementary/complete determinant duality,
  and the second-kind Stirling moment determinant agree throughout this box.
- The chain/subset bijection agrees on 11,956 enumerated chain instances
  across these triples; this count includes repetition over k.

The final output line is `ALL_CHECKS_OK`; the program exits 0. These finite
checks support transcription and the proposed mathematical route. They do
not establish the unbounded conjecture or claim a Lean proof.
