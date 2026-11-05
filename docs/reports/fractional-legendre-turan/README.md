# Fractional Legendre Turán comparison

`check.py` checks the strict inequality in Dmitry Karp's SF21 Conjecture 4:

\[
F_\mu(y)^2>F_{\mu-\delta}(y)F_{\mu+\delta}(y),\qquad
F_\nu(y)={}_2F_1(-\nu,1+\nu;1;y),
\quad \mu>0,\quad 0<\delta,y<1.
\]

Source: [NIST-hosted SF21 slides](https://math.nist.gov/~DLozier/SF21/SF21slides/Karp.pdf),
physical PDF pages 52–54, conference 6–8 April 2011. The source's preceding
integer-degree result supplies the argument interval. No positivity condition
on the values of the functions is included in the conjecture.

Run with Python 3 and SymPy (checked with SymPy 1.14.0):

```sh
python3 docs/reports/fractional-legendre-turan/check.py
```

`output.txt` is the output of this command. All calculations use integers,
`fractions.Fraction`, or exact SymPy expressions. There is no floating-point
special-function evaluation.

The fidelity checks recover the source's integer-degree identity
\(T_{1,1}(y)=2y(1-y)>0\), and verify 936 rational grid values with integer
central degrees 1 through 12. Eighteen terminating-series values independently
check the enclosure implementation.

The refutation checks enclose three noninteger-degree values of the original
series, including

\[
-\frac{1256997}{100000000}\le
T_{3/2,127/128}(255/256)
\le-\frac{314249}{25000000}<0.
\]

All three degrees are positive. The function values themselves have mixed
signs: the central value is positive and both neighboring values are negative.

The enclosure algorithm starts with the constant term 1 and uses the exact
term ratio

\[
\frac{u_{j+1}}{u_j}
=y\frac{(j-\nu)(j+\nu+1)}{(j+1)^2}
=y\left(1-\frac1{j+1}-\frac{\nu(\nu+1)}{(j+1)^2}\right).
\]

Each update uses outward integer floor/ceiling at a fixed dyadic scale.
The interval sums therefore enclose the finite partial sum, including its
negative early terms. Once \(N>\nu>0\), all subsequent ratios are positive
and smaller than \(y\). The remaining absolute tail is at most
\(\lvert u_N\rvert y/(1-y)\). Interval multiplication and squaring then enclose
the Turán difference. The printed coarse rational intervals are also rounded
outward. These are rigorous real-series enclosures, independently of any
endpoint continuation formula.

The symbolic checks verify the neighboring-digamma recurrence algebra,
the trigonometric product identity, and the negative linear coefficient
\(-2a_\mu^2/(\mu(\mu+1))\) in the proposed asymptotic family. They do not
certify DLMF 15.8.10, its uniform remainder, or a universal asymptotic theorem.
Those are analytic proof obligations in the preregistration. The finite
counterexamples already contradict the source assertion; a Lean proof is
outside this experiment.

Analytic ingredient: [DLMF 15.8.10](https://dlmf.nist.gov/15.8.E10),
with integer parameter zero. The checker is original source code; no article
text or third-party implementation is copied into this report.
