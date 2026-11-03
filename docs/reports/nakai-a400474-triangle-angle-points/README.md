# Nakai A400474: exact triangle-angle count checks

Source: Takuto Nakai, [OEIS A400474](https://oeis.org/A400474), revision 12
(2026-11-02 14:26:30, OEIS displayed timestamp), submitted 2026-10-25.
The entry calls the six-residue formula a conjecture. The copied Python function
and Karan Sharma b-file are source material under the
[OEIS End-User License Agreement](https://oeis.org/LICENSE).
The b-file snapshot is `oeis-b400474.txt`.

Run with Python 3 and no third-party packages:

```sh
python3 check.py
```

`output.txt` is the captured output; the final line is `ALL_CHECKS_OK`.
All calculations use integers or `fractions.Fraction`.
The source program is reproduced for n = 3 through 30, all 998 b-file entries
are compared with the conjectured formula, and an independent enumeration of
distinct positive triples checks the full floor-threshold optimization for
n = 3 through 199. Direct point-triple filters are compared for n <= 22.
The n = 180 examples are checked through the source's discrete reduction.

For ordered distinct positive triples u + v + w = n, the exact number removed
by u >= n-r, for 0 <= r < n, is

    H(n,r) = floor((r-1)^2/2) - 2 max(0,r-floor(n/2))
             + 2 [3 divides n and 3r >= 2n].

The simpler F(r) = floor((r-1)^2/2) applies when r <= floor(n/2).
For example H(6,5) = 6 whereas F(5) = 8. The checks cover this distinction,
the donor/receiver inequality used to reduce large thresholds, and the
balanced minima in all six residue classes.

These are finite exact checks of the discrete model. They neither prove the
geometric bijection from Euclidean interior points nor prove the conjecture
for all n. The proposed settlement uses a general geometric argument and
an unbounded exchange/convexity proof; enumeration is supporting evidence.
