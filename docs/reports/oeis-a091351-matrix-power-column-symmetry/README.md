# OEIS A091351: matrix-power column symmetry

Paul D. Hanna's final FORMULA conjecture in [OEIS A091351](https://oeis.org/A091351), published revision 37 (2026-07-17), asserts

\[
(T^p)_{n+q-1,q-1}=(T^q)_{n+p-1,p-1}
\qquad(n\ge0,\ p,q\ge1).
\]

The matrix is lower triangular, indexed from zero, and defined by requiring T(n,k) to equal row sum n-k of T^k.

Run with Python 3, without third-party packages or network access:

```sh
python3 check.py
```

The stored output is [check.out](check.out). All 66 published DATA terms and all 63 displayed entries of powers 2, 3 and 4 agree. Three independent constructions (the NAME row-sum definition, the first FORMULA scalar recurrence and the quadratic-shift recurrence) agree on a 35 by 35 truncation. All 630 defining row-sum equations and 3,962 symmetry cases with 1 <= p,q <= 12 hold. The example n=4, p=2, q=3 has common value 845.

The program also checks the proposed finite algebraic route for sizes 1, 2, 3, 5, 8 and 12: nilpotent polynomial inversion, conjugation, the factorization of S^m T^m, the first-vector identity and commutation. These checks use exact Python integers throughout. They provide finite evidence and regression cases, without proving the assertion for unbounded indices or establishing a Lean proof.

The DATA and POWER_ROWS constants reproduce material from Paul D. Hanna and the OEIS Foundation, Inc., under the [OEIS CC BY-SA 4.0 terms](https://oeis.org/LICENSE). The source entry and its indices, including both shifts in the conjecture, are retained in the program's citations and checks.
