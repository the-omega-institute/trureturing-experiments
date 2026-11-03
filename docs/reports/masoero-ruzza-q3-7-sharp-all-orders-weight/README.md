# Masoero–Ruzza Question 3.7 exact check

This report checks Question 3.7 of Masoero–Ruzza, arXiv:2605.24563v3. The Python 3 program constructs the source's Laguerre Wronskian from the exact monic generalized Laguerre coefficients, extracts the monic polynomial \(\mathcal P_\lambda^{(\beta)}\), and computes root power sums by Newton identities.

All arithmetic is exact: SymPy rational polynomial domains and interpolation at rational half-integer values; no floating point is used. It reproduces the source's \(S_1\) identity and its displayed conjectural \(S_2\) formula for every partition of size at most four at 17 rational beta values. It then interpolates \(S_k\) for \(k=1,2,3,4\) and checks degree at most \(2k+2\), with independent exact evaluations at three fresh beta values, for the same 15 partitions.

Run with:

```text
python3 check.py
```

The final line is `ALL_CHECKS_OK` on success.
