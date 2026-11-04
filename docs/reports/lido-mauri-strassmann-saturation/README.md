# Lido–Mauri saturation stopping criterion

Source: Guido Maria Lido and Luca Mauri, *A multivariate Strassmann theorem*, arXiv:2605.03683v1, Question 3.2.4 and Example 3.2.3. These are independent exact checks supporting a preregistration; they are not a Lean proof.

Run with Python 3 and SymPy:

```sh
python3 docs/reports/lido-mauri-strassmann-saturation/check.py
```

The source-fidelity check reproduces the colon certificate, the three lifted syzygies, and the reduction in Example 3.2.3, treating its unknown series as formal symbols in universal polynomial identities. The counterexample check verifies five polynomial identities over the integers, the finite quotient with basis `1,z,w`, and the square-zero algebra obstruction over the rationals. The identities and the nonzero coefficient −2 give the same argument over every coefficient field of characteristic different from 2. The generic point `(y,z,w)=(0,0,1)` verifies nondegeneracy.

The saved output is `check.out`; its last line is `ALL_CHECKS_OK`. All calculations use exact integers, rational numbers or symbolic polynomials. The mathematical argument must still justify the restricted-series evaluation map, reduction of the colon ideal, and the Krull-dimension conclusion; this program does not certify those bridges in Lean.
