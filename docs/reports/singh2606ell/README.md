# Singh unital-qubit weighted ellipsoid check

This directory checks the conjecture after Proposition 4.2, equations
(4.11)–(4.12), in Satvik Singh, *Gaussian mean width strong converse bound on
the classical identification capacity of quantum channels*,
arXiv:2606.05032v1 (3 June 2026).

Run with Python 3.12 or newer and SymPy 1.14:

```sh
python3 docs/reports/singh2606ell/check.py
```

`check.out` is the complete successful output. The program exits zero and ends
with `ALL_CHECKS_OK`. It reproduces the mixed-state Gram normalization and
second-moment bound in the source. It checks the determinant-scaled singular
value prefix and positive power-sum inequalities proposed for the settlement.
The script is independently written from matrix formulas in the paper and the
proposed mathematical argument; no third-party program source is copied.

The complete small grid consists of all nonnegative integer Pauli-mixture
4-tuples with total 1 through 4, at each of three reference-state Bloch axes and
three polarizations. The 1,200 random cases (seed 439) use Pauli mixtures with
independent quaternion rotations of input and output and arbitrary rational
reference-state orientations. Ranks one through four and unitary channels are
included. There are 1,821 cases, 7,284 prefix comparisons, and 16,389 power-sum
comparisons at exponents 1/10, 3/10, 7/10, 1, 3/2, 2, 3, 5 and 10.

All evidence uses exact rational matrices and certified rational root bounds.
For a positive Hermitian matrix K, sigma = K²/Tr(K²) gives a rational weighted
superoperator even when sigma's inverse square root is irrational. Characteristic
roots are isolated using exact polynomial intervals. Rational powers use integer
nth-root bounds whose defining inequalities are checked exactly. Six spectral
prefix comparisons that straddle equality are certified by exact compound-matrix
determinants. No floating-point value is accepted as evidence.

The tested algebraic identities include J W J = c² W⁻¹ and, in the invertible
cases, G C⁻¹ G = J C J. The finite checks detect errors; they do not prove a
universal conjecture or certify a Lean proof. The stored run used Python 3.14.8
and SymPy 1.14.0.
