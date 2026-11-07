# Interior unistochasticity check

`check.py` is an independent exact-arithmetic check for the Rajchel–Gąsiorowski–Życzkowski 2018 Question 3 preregistration. It reproduces a source triangle result at order 4 with an exact orthogonal lift, verifies the cubic Sidon support in `F_32 = F_2[X]/(X^5+X^2+1)`, checks the trace-word dimension inequality used by the proposed family refutation, and verifies the skew-Hadamard recursion through order 128. The script uses only integers, `fractions.Fraction`, and SymPy exact algebra; it uses no floating point.

Run:

```text
python3 check.py
```

The final line must be `ALL_CHECKS_OK`. These checks do not enumerate unitary lifts and do not constitute a Lean or kernel verification of the proposed theorem.
