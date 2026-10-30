# Sun A346934: composite-modulus clause

Exact finite checks for the second sentence of Conjecture 1 in
[OEIS A346934, revision 29](https://oeis.org/A346934): the permanent of
`[i-j]` of order `2n` is divisible by `2n+1` when that modulus is composite.
The sign sentence and Conjecture 2 are excluded.

Run with Python 3, without third-party dependencies:

```sh
python3 docs/reports/sun-a346934-composite-modulus/check.py
```

`output.txt` contains the complete output. The program reproduces the source's
example and eleven sequence terms and the prime-modulus residue stated in
Sun, *Arithmetic properties of some permanents*, arXiv:2108.07723v7,
Theorem 1.2(ii). Independent permutation summation and subset dynamic programming
check the rank-two expansion, including signed and degenerate inputs. It checks
every factorial split at every odd composite modulus through 999, the nominated
conjecture for `1 <= n <= 100`, and arbitrary signed additive matrices at moduli
9 and 15. Random inputs use the fixed seed 346934.

All evidence uses integer arithmetic. These finite checks support a proposed
proof route; they do not establish the universal conjecture or any Lean theorem.
The generic expansion is a known rank-two permanent identity, not claimed as new.
