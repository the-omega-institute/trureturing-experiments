# Chemical antiregular Conjecture 10

Source: M. Knor, R. Škrekovski, S. Filipovski and D. Dimitrov,
*Extremizing antiregular graphs by modifying total σ-irregularity*,
arXiv:2411.01530v1, Conjecture 10; Applied Mathematics and Computation
490 (2025), 129199, DOI 10.1016/j.amc.2024.129199.

Run with Python 3, without third-party packages:

```sh
python3 docs/reports/knor-skrekovski-antiregular-conjecture-ten/check.py
```

The program uses only integers. It bounds each positive nth root of 2 or 3
between consecutive multiples of 1/2^96 and verifies both endpoints by
integer exponentiation. Signed linear combinations use outward bounds.
Every nonmaximizing vector has a certified strictly positive gap; tied
maximizers have identical coefficient triples.

The source-fidelity check reproduces Theorem 9 at n=8 with exponent 1/64.
Its hypothesis is certified by the integer inequality
3 × 23^64 < 24^64. The settlement check enumerates all 407,363
nonnegative degree-count vectors with even degree sum for 7 ≤ n ≤ 64,
including nongraphical vectors, and finds precisely the Conjecture 10
maximizers for exponent 1/n. It checks the quadratic-coordinate identity
coefficient by coefficient, direct pair sums on five connected bases,
and 25 degree-preserving insertions per base.

`output.txt` contains the complete output, including explicit base edges
and rational lower bounds on each strict optimization gap. Its final line
is `ALL_CHECKS_OK`. These finite checks support a planned argument for
every n ≥ 7; they do not prove the unbounded conjecture.

The universal argument consists of a positive quadratic penalty on a
congruence lattice, a uniform bound on the displacement of its scalar
maximum, and a connected insertion adding one vertex of each degree
1, 2, 3 and 4. Positive finite enumeration is only a check of this argument.
