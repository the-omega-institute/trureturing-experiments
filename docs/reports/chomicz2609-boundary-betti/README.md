# Secondary boundary Betti elements

`check.py` checks the gcd criterion answering the unnumbered question after
Lemma 3.13 in §3.2 of K. Chomicz, arXiv:2609.04000v1 (3 September 2026),
*The type and cardinality of minimal presentations of numerical semigroups
with embedding dimension four*.

Run with Python 3.12 or newer:

```sh
python3 docs/reports/chomicz2609-boundary-betti/check.py
```

All calculations use exact integers. `check.out` is the complete output.

The fidelity check reproduces Example 2.9: the primary critical coefficients
are (17,7,15,4); the secondary example (145,203,74,111) has coefficients
(7,5,3,2) and exactly the four factorization fibers displayed in the paper.

The exhaustive check considers every increasing tuple of four integers in
[2,40], filters by gcd one and minimal generation, and tests all three
partitions into critical pairs. Lemma 2.2 of the source ensures every
secondary semigroup occurs in this search. Each retained tuple is also
checked directly against the primary/secondary definition using its
critical fibers. All eight source-compatible generator orderings are
checked, for both i=2 and i=3.

Among 28,106 minimally generated gcd-one tuples, 85 are secondary. Their
680 compatible orderings give 1,360 boundary cases: 196 are Betti and 1,164
have connected factorization graphs. In every case, with p=gcd(d0,d1),
b_ii*d_i is Betti exactly when b_ii=p/gcd(p,d_i). In the affirmative cases,
the two components are precisely the factorizations supported on {0,1}
and on {2,3}. A reversed-loop factorization enumerator and a separate
union-find graph calculation cross-check every boundary case.

The source example has b22=16 and b33=11, with gcd threshold 29 for each;
both graphs are connected. A positive example is (16,24,10,15), i=2:
b22=4 equals its threshold, and degree 40 has two components. A negative
example is (10,15,14,21), i=2: b22=4 is below threshold 5, and degree 56 has
one component.

These are finite checks supporting a proposed general argument. They do
not establish the unbounded claim or constitute Lean proof evidence.
