# Bala's A007820 square-divisibility conjecture

Source: [OEIS A007820](https://oeis.org/A007820), revision 89. The entry defines
`a(n) = S(2n,n)` and records Peter Bala's conjecture that `a(2n+1)` is divisible
by `(2n+1)^2`. It also records Vladimir Kruchinin's finite reflection identity.

Run with Python 3 and the standard library:

```sh
python3 check.py
```

The program computes Stirling numbers by their integer recurrence. It reproduces
all 18 displayed source DATA values, checks inclusion-exclusion for 676 pairs,
and checks the source reflection identity for every index from 0 through 201.
It checks square divisibility and the first-order-term extraction for the 101
positive odd indices through 201. The even index 2 provides a boundary check.
All arithmetic is exact; no floating point or external package is used.

`output.txt` is the program's standard output. Its final line is `ALL_CHECKS_OK`.
These finite checks support source fidelity and the proposed algebraic route;
they do not prove the universally quantified conjecture.
