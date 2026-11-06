# OEIS A397594: quadrupled-row dyadic parity

`check.py` uses only exact Python integers and `fractions.Fraction`. It recursively solves

\[
[x^n]\frac{A(x)^{n+1}}{1-(n+1)x}=\bigl(4(n+1)\bigr)^n
\]

for the coefficient of degree `n`, then checks the equivalent rescaled row sum from the OEIS definition. It reproduces the first 25 values of [A397594's b-file](https://oeis.org/A397594/b397594.txt), verifies both source identities on finite ranges, and checks the proposed odd-support set through degree 100. The finite checks are source-fidelity and counterexample screening only; they do not prove the unbounded conjecture.

Run:

```text
python3 check.py
```

The successful output is stored in `output.txt`.
