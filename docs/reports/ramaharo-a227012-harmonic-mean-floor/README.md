# Ramaharo's cubic-block harmonic-mean floor conjecture

The target is the final FORMULA line of [OEIS A227012](https://oeis.org/A227012),
revision #45 (2024-11-03), attributed to Franck Maminirina Ramaharo (2018-04-12).
For every integer n ≥ 1 it predicts

\[
\left\lfloor\frac{n^3-(n-1)^3}
{\sum_{k=(n-1)^3+1}^{n^3}1/k}\right\rfloor
=\frac{2n^3-3n^2+n+2}{2}
+\left\lfloor\frac{\max(0,n-3)}4\right\rfloor.
\]

Run with Python ≥ 3.12 and SymPy ≥ 1.14:

```sh
python3 check.py > output.txt
```

The program reproduces all 42 published data values from exact rational
reciprocal sums and checks the formula for n=1,...,64, including the exact
initial harmonic means 1 and 1960/481. It checks the central moments through
degree six symbolically, the symmetric reciprocal identity and finite geometric
remainder, and both uniform rational certificates. Their numerators have
strictly positive coefficients in z=n−3. It also checks the four residue classes
used for rounding. No floating point arithmetic is used.

With p=(n−1)³, q=n³, N=q−p, c=(p+q+1)/2, h=(N−1)/2 and
P=n³−3n²/2+3n/4+3/8, the certificate is designed to establish
P−1/8 < c/B ≤ N/H ≤ c/A < P+1/8 for n≥3. Here
A=1+μ₂/c²+μ₄/c⁴ and B=A+μ₆/(c⁶(1−h²/c²)); the μ's are
the exact central moments printed by the program. The unbounded inequalities
still require their mathematical argument and Lean verification; this check is
not a Lean proof. Kimberling's broader polynomial-block recurrence conjecture
is outside the target.
