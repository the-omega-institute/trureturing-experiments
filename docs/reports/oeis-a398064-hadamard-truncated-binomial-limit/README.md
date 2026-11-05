# OEIS A398064: truncated binomial Hadamard limit

`check.py` uses exact rational arithmetic and SymPy algebra. Run with Python 3.10+ and SymPy 1.14 (`python3 -m pip install sympy==1.14.0`):

```sh
python3 check.py
```

The check reproduces the equality of the source's digamma and radical/logarithm expressions for H(1/4). It verifies symbolic limits at every integer x from -8 through 8, exact finite-sum identities, rational bounds for the x=-1 and x=-2 limits, and finite coefficient identities at rational and complex parameters. It uses no binary floating point or decimal approximation. `check-output.txt` contains the resulting output.

For negative integers the exact polynomial summation leaves the harmonic block H_(2n)-H_(n-1). Its limit log(2) is the standard integral/Riemann-sum fact; the program separately checks that the block's multiplier and the remaining rational part have the asserted finite symbolic limits. The rational log(2) enclosure uses its positive atanh(1/3) series and an explicit geometric remainder.

These checks support a preregistration. They do not prove the limit for arbitrary complex x, certify 105 decimal digits, or supply a Lean proof. The independent analytic argument and the binding statement belong in the preregistration draft.

Source: David Ulgenes, [OEIS A398064](https://oeis.org/A398064), revision 39; the FORMULA line beginning “Conjecture: H(x+1)”. OEIS source material is subject to the [OEIS End-User License Agreement](https://oeis.org/LICENSE). The program implements the mathematical formulas independently and includes no copied source code.
