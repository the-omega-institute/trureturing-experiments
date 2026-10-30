# Hanna A381355: exact preregistration checks

Paul D. Hanna's [OEIS A381355](https://oeis.org/A381355) conjectures that
its n-th term is divisible by the n-th prime for n > 1. Its source series
is [A381353](https://oeis.org/A381353): F(0) = [x]F = 1 and
[x^n]F(x)^prime(n) = 0 for n > 1. The target term is [x^(n-1)]F'/F.
A381355 revision 13 and A381353 revision 19 supply the quoted data.

Run with Python 3, using only its standard library:

```sh
python3 docs/reports/hanna-a381355/check.py
```

The rational triangular construction is computed before integrality is tested.
Truncated convolution and binary powering reproduce both DATA prefixes and
three displayed coefficient-table rows. The source is checked through degree
64. Eighteen combinations of integer slopes and prime schedules are checked
through degree 24. An independent inverse-series calculation and the
logarithmic-derivative recurrence agree. Diagonal zero equations, low-degree
prime residues, the derivative identity, rational coefficient perturbations,
and prime divisibility all pass. The excluded n = 1 case fails divisibility,
as the source permits. Full output is in `check-output.txt`.

These finite checks support fidelity and the proposed proof route. They do
not establish the universal conjecture, priority, or a Lean theorem.

The independently written check implementation is provided under the MIT
license. The short OEIS DATA excerpts in the program are attributed to
Paul D. Hanna and OEIS and retain [CC BY-SA 4.0](https://creativecommons.org/licenses/by-sa/4.0/)
terms. No source paper, downloaded snapshot, or research transcript is included.
