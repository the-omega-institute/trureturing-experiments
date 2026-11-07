# Kang–Lee–Lim single-lamination log-concavity check

This report records exact arithmetic checks for Conjecture 1 of Kang, Lee and Lim,
*Unimodality and Cluster Algebras from Surfaces*, arXiv:2508.04396v3, Section 5.1.
`check.py` uses SymPy symbolic polynomials and integer matrix mutation. It verifies
the source's octagon example (`2q^3 + 4q^2 + 2q + 1`) and the nominated torus
witness, whose two mutations specialize to `(1+q)^2 + q^3`. It then checks the
family `(1+q)^2 + q^m` for `m = 3` and `m = 4`; at the failing index `m-1`, the exact inequality is
`0 = a_(m-1)^2 < a_(m-2) a_m = 1`.

Run with Python 3:

```text
python3 check.py > output.txt
```

The check uses no floating-point values. The two nominated values are a reproducibility
check; the proposed mathematical settlement will require a formal geometric
certificate for the genus-one single lamination and the two mutations.
