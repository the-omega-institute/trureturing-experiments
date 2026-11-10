# Affine panmagic span checks

Run with Python 3.12 or later from the repository root:

```sh
python3 -u docs/reports/koshkin-lee-2025-affine-panmagic-span-complete-classification/check.py
```

`check.out` is the program's complete output. A successful run exits zero and
ends with `ALL_CHECKS_OK`. The self-authored program uses only the Python
standard library and exact integer arithmetic.

For a Pólya order n > 1, let A_n be the real span of the permutation matrices
of x ↦ ax+b modulo n, where a, a−1 and a+1 are units. Let P_n be the real
space whose row, column and both toroidal diagonal-family sums all have a
single common value. The proposed uniform classification is A_n = P_n
exactly at prime orders, with

    dim A_n = product over p^e || n of
              (1 + (p−1)(p−3) * sum(p^(2*j), j=0,...,e−1)).

The source is Koshkin–Lee, *Panmagic permutations and N-ary groups*, PUMP
Journal 8 (2025), 195–212, DOI 10.46787/pump.v8i.4862, arXiv:2606.22221v1,
Section 8's unnumbered concluding affine spanning question.

The checks reproduce the source's order-five statement: all ten panmagic
permutations are affine, their matrices have rank nine, and summing the five
intercept translates at either slope yields the all-ones matrix. They also
check the three-unit characterization for every affine permutation at orders
2 through 60.

For all 49 Pólya orders between 5 and 150, the program enumerates the exact
Fourier support, checks its CRT decomposition, checks the prime-power product
formula, and compares it with the full panmagic support. Separate integer row
elimination computes rational ranks directly from the permutation matrices
and from the literal equal-line-sum constraints:

| n | affine rank | panmagic dimension |
|---|---:|---:|
| 5 | 9 | 9 |
| 7 | 25 | 25 |
| 11 | 81 | 81 |
| 13 | 121 | 121 |
| 25 | 209 | 529 |
| 35 | 225 | 1089 |

At each of the 16 composite Pólya orders in that interval, an integer
Ramanujan witness R(x,y)=c_n(px+y), using the least prime divisor p, has zero
sums in every line family, zero inner product with every affine generator,
R(0,0)=φ(n)>0, and squared norm n²φ(n). Its positive perturbation of the
normalized constant matrix is checked after clearing the denominator.

Integer elimination uses only nonzero row scaling, subtraction and removal
of a common gcd. It computes rational rank, which equals real rank for these
integer matrices. The witness uses the exact divisor formula for Ramanujan
sums; no floating-point trigonometry or modular-rank inference is used.

These are finite consistency checks. The result for all Pólya orders requires
the uniform Fourier-support, prime-power and CRT arguments, and the uniform
missing-frequency witness; this program supplies neither a Lean proof nor
an exhaustive literature certification.
