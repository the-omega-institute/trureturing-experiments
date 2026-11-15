# Dagar Conjecture 1.1: exact Mackey checks

Source: Prem Dagar, *Parabolically induced representations of GL_n(O_ell)*,
arXiv:2610.15762v1, Conjecture 1.1 and Section 3.1.

Run with Python 3.12 or newer, using only its standard library:

```sh
python3 docs/reports/dagar2610c11/check.py
```

`check.py` covers O = Z_p, p in {2,3}, n in {2,3}, ell = 2,
and every multiplicative-character tuple (264 in total). Matrices have integer
entries modulo p^2. The unit groups are cyclic, of orders 2 and 6; character
values are represented exactly by root-of-unity exponents. No floating point
arithmetic is used.

Upper triangular column elimination canonically represents right B-cosets.
Generating G by elementary matrices and diagonal unit generators enumerates G/B.
The number of flags is checked against the exact group and Borel orders, and
for n=3 against the source's dimension formula
p^3 (p+1) (p^2+p+1). B-orbits on these flags enumerate B\G/B. Schreier generators
of each orbit stabilizer supply all diagonal-character compatibility constraints.
The Mackey inner product is the number of compatible double cosets. For n=2
and for p=2, n=3, direct enumeration of B intersection conditions independently
checks every computed inner product.

`output.txt` lists all tuples and their exact inner products. The respective
numbers of double cosets are 3, 18, 3, 18. Every inner product is 1 exactly when
the character exponents are pairwise different modulo p, equivalently when their
restrictions to 1+p Z/(p^2) are pairwise distinct. There are 74 irreducible cases.
The rank-three check reproduces the source's irreducibility theorem for these
parameters, including residue characteristic 2. The output ends in
`ALL_CHECKS_OK`. These finite checks support the general proposed argument;
they do not prove the all-rank, all-length conjecture.
