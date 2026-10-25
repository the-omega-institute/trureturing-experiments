# Provvedi's prime-period conjecture for OEIS A249938

Source: [OEIS A249938](https://oeis.org/A249938), revision 76. The target is
Federico Provvedi's second conjecture, restricted to its entire prime-modulus
clause and the positive-index convention of its example. The separate
288-divisibility assertion and composite moduli are outside the target.

Run `python3 docs/reports/provvedi-a249938/check.py` from the repository root.
Python 3.8 or later suffices; there are no third-party dependencies. The program
uses integer arithmetic exclusively and raises an exception on any failed check.
`output.txt` is its complete output, ending in `ALL_CHECKS_OK`.

The check constructs Stirling numbers by their integer recurrence and forms
`F_N = sum(k! S(N,k))`. An independent ordered-partition recurrence agrees through
`F_80`. The initial values and the source's eleven residues modulo 23 are reproduced.
For all 46 primes from 2 through 199, it checks the predicted period, excludes
all smaller cyclic shifts, and excludes smaller shifts in four finite tail windows.
It also checks the full finite-field exponential formula, opposite-residue pairing,
zero-weight count, square-frequency support, and nonzero Vandermonde determinant.
The integer sequence is computed through positive index 413. The index-zero
obstruction for every odd prime tested is checked separately.

These are finite checks. They do not prove the all-prime assertion, any infinite
tail assertion, or prior-art novelty. No Lean proof is included.
