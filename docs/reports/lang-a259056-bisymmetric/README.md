# Lang's bisymmetric determinant conjectures

`check.py` uses Python 3's standard library and exact integers and rational
numbers. Run it from this repository's root:

```sh
python3 docs/reports/lang-a259056-bisymmetric/check.py
```

The stored result is `output.txt`. The program constructs each matrix from
Lang's literal consecutive wedge and its symmetry orbit, reproduces the
four displayed matrices and the twenty determinant terms of
[OEIS A259056 revision 7](https://oeis.org/A259056/internal), and reproduces
the fifteen quotient terms of
[OEIS A259057 revision 7](https://oeis.org/A259057/internal).

For every dimension 1 through 32 it compares fraction-free Bareiss elimination
with independent rational Gaussian elimination, verifies an explicit reflection
basis and its inverse, checks both reflection blocks and their triangular
min-kernel factorizations, and checks the odd factorial identity and even double
factorial identity. Every checked even determinant is divisible by three. The
empty min-kernel case is included. The comparison sequence is
[OEIS A034976 revision 30](https://oeis.org/A034976/internal).

These are finite checks supporting a preregistration draft. They provide no
all-dimension proof, and no Lean proof has been written or checked here.
