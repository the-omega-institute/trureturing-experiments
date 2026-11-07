# Prescribed 14325 packing counts modulo 10

Source: Natasha Blitvić, Mohammed Slim Kammoun and Einar Steingrímsson,
*Permutations, Pattern Covers, Posets and Positivity: Part I: Patterns and
Posets*, [arXiv:2610.29846v1](https://arxiv.org/abs/2610.29846v1), Section 7
item 11; Definitions 2.1–2.3, 2.5–2.6 and Remark 2.7. The source asserts an
unexplained residue pattern through 23 occurrences. Only the modulo-10
observation is considered here.

For positive n, A_n counts permutations of length 3n+2 whose consecutive
five-entry windows at starts 1,4,...,3n-2 have relative order 14325.
Additional occurrences are permitted. The proposed assertion is A_n mod 10
= 1 if n+1 is a power of two, and 6 otherwise.

Run with Python 3; there are no third-party dependencies:

```sh
python3 check.py
```

`output.txt` is the program's exact output. A failing condition raises an
exception and a nonzero exit; success ends with `ALL_CHECKS_OK`. Arithmetic
uses integers only.

The check reproduces Remark 2.7's extra-occurrence example and all 16 exact
counts displayed in item 11. It independently enumerates actual permutations
for n=1,2, then counts ideals of the poset obtained directly from **all**
window comparisons for n=1,...,10. These original counts agree with an
independently built grid and the weighted recurrence
B_1(k)=1, B_n(k)=sum_{i=1}^k i B_{n-1}(i+2), A_n=B_n(3).
This recurrence is credited to Michael Wallner in OEIS A274644; the grid
decomposition is already in Ran Pan's *Algorithmic Solution to Problem 1*
(2016). Neither ingredient is claimed as new.

The program also compares 30 generalized grids (n=1,...,5, first-rung
internal length v=0,...,5) with B_n(v+2). For n=1,...,256 it checks the
modulo-5 residue vectors for k=0,...,65, the ballot-triangle parity for
r=0,...,20, the Catalan identity and parity, and the proposed modulo-10
formula. The exceptional n in that range are 1,3,7,15,31,63,127,255.

These are finite supporting checks. They do not prove the source-to-grid
bijection, the recurrence for arbitrary n, or the universal congruence.
The preregistration's later Lean probe must establish those bridges without
assuming a recurrence in place of the original permutation count. No Lean
proof of the assertion or its negation is included.
