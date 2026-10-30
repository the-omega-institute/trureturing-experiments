# Permuted rotor loop without the star condition

Source: Nik Jenzer, *Star-Ordered Rotor Loops in Clifford Algebras*,
arXiv:2610.21593v1, §5.4, Theorem 5.5 and the following scope remark
(Remark 5.7 in the TeX, rendered Remark 5.5 on page 19).

Run with Python 3, using only the standard library:

```sh
python3 docs/reports/permuted-rotor-loop/check.py
```

`output.txt` is the program's output. The program exits 0 and ends with
`ALL_CHECKS_OK`. All arithmetic uses integers and exact rational numbers.
Sparse polynomial multiplication retains independent variables `c_i,s_i`,
so the checks include zero coefficients and both signs. For the source's
angles, substitute `c_i=cos(theta_i/2)`, `s_i=sin(theta_i/2)` and
`c_i^2=1-s_i^2`.

- Fidelity: reproduces Remark 4.4's triangle in `Cl_{0,3}`:
  `B1*B2*B3=-1` and `scalar(P3)=c1*c2*c3-s1*s2*s3`.
- Checks all six return orders for that triangle and for an exact rational
  rotation of its bivectors. Both families violate the star condition.
- Checks every permutation for `1<=k<=6` (873 orders) by direct symbolic
  multiplication in the abstract anticommuting presentation. It separately
  checks the odd full-support coefficient, which a tracial functional on
  a representation can retain, and finds it exactly zero. Explicit finite
  conjugation witnesses exclude every other nonempty trace component.
- Checks all 142 canonical bivector subset families of square `-1` with
  pairwise anticommutation in every nondegenerate signature in dimensions
  2 through 4, with all 326 return orders.

These finite exact checks support the expected settlement **Proved**.
They do not certify the unbounded claim. No Lean proof is included.
