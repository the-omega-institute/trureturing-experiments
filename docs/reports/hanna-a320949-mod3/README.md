# Hanna A320949: exact finite checks modulo three

The source is Paul D. Hanna's [OEIS A320949](https://oeis.org/A320949), revision 12, whose FORMULA line conjectures `a(n) = 0 (mod 3) for n > 0`. Its defining equation is

```text
Sum_{j>=0} (1+x)^(j^2) / (1+A(x))^(j+1) = 1,  A(0)=1.
```

Run with Python 3 and its standard library:

```sh
python3 docs/reports/hanna-a320949-mod3/check.py
```

`output.txt` contains the output. The program uses integers and `fractions.Fraction`; its evidence contains no floating point arithmetic. It independently computes the finite transform of each coefficient in two ways: Newton differences of `binom(j*j,k)`, and exact polynomial moments of the geometric weights `2^(-j-1)`. All 325 transform coefficients with `k+r <= 24` agree. The two resulting recursions reproduce `a(0)..a(24)` from the [source b-file](https://oeis.org/A320949/b320949.txt), and every nonconstant term is divisible by three.

The program also tests the triangular coefficient dependence on an unrelated integral prefix, and the rational-function reciprocity identity and its mod-three consequence for `binom(t*t,k)`, `k=0..60`, at four rational arguments. These finite checks corroborate a proposed proof route; they establish neither the unbounded conjecture nor a Lean proof. The last line is `ALL_CHECKS_OK` on success, and a failed assertion exits nonzero.

The source values are quoted with attribution under the [OEIS End-User License Agreement](https://oeis.org/LICENSE). The check program follows this repository's code license.
