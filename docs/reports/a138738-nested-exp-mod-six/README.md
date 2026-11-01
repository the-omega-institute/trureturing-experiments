# Exact checks for Hanna's nested-exponential congruence

Source: Paul D. Hanna, [OEIS A138738](https://oeis.org/A138738), revision 6.
The entry defines the EGF by nested exponentials with successive powers
x, x², x³, … and conjectures that a(n) is 1 modulo 6 at every index.
The 22 reference terms and displayed EGF coefficients in `check.py` come
from that entry; see the [OEIS license](https://oeis.org/LICENSE).

Run with Python 3 (standard library only):

```sh
python3 docs/reports/a138738-nested-exp-mod-six/check.py
```

`check.out` is the captured output. The verifier uses integers and
`fractions.Fraction` throughout. It reproduces all 22 displayed terms and
the displayed EGF through degree 5, checks integrality and the congruence
through index 40, compares truncation depths, and independently computes
exponentials by a direct sum of powers through degree 24. It also compares
the depth-multiplicity formula through index 20 (490 weighted-partition
terms), checks its factorial factors, and checks the EGF differential
recurrence through index 39. All checks finish with `ALL_CHECKS_OK`.

These bounded computations support source fidelity and the proposed
structural argument. They do not establish the unbounded conjecture and
contain no Lean proof or refutation.
