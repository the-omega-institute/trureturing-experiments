# HMR universal nonuniform-reflection belief

Exact arithmetic checks for the universal belief immediately after Conjecture
3.13 in U. Haagerup, M. Musat and M. B. Ruskai, *Extreme Points and
Factorizability for New Classes of Unital Quantum Channels*,
[arXiv:2006.03414v4](https://arxiv.org/abs/2006.03414v4), Section 3.3.2,
using the matrices of equation (10). The existential Conjecture 3.13 is a
different statement and is not settled by these checks.

Run with Python 3 and SymPy:

```sh
python3 check.py
```

The program uses integer and rational arithmetic and exact symbolic rational
functions. `check.out` is its complete successful output. It checks:

- Source fidelity: the uniform-reflection product family has full rank for
  dimensions 3, 4 and 5 at parameter zero, as in Theorem 3.11(a), and has the
  normalization stated before equation (11). A nonuniform dimension-three
  sample also has full product rank at the exceptional parameter.
- The rational family
  (y(z)=(0,2,-6z,9z^2-1,6z,2)/(3(1+3z^2))) has unit norm,
  orthonormal cyclic translates and nonuniform nonzero last five coordinates
  for rational (z\notin\{0,\pm1/3\}).
- For its source matrices, the exact symbolic identities are
  (sum_m A_m=(t-3)I), (sum_m A_m^2=(5+t^2)I), and
  (sum_{m\ne n}A_mA_n=(4-6t)I). They give a nonzero dependence at
  (t=-1/5).
- For 36 rational parameters, including (z=1), it computes the source's
  matrices literally by cyclic conjugation, checks every coefficient in the
  dependence, and checks that the six matrices are independent and invertible.
  At (z=1), the ordered product rank is 22 of 36.

The symbolic identities are algebraic evidence for the whole family; the
36 rational samples provide a separate finite check. This report contains
no Lean proof and makes no literature-priority claim.
