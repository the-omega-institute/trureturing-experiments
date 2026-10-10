# Mixed-sieve cumulative response

The program evaluates the fixed weights used in the proposed local
Robin negative-part estimate. It uses the arbitrary parameter samples

\[
A=m^8,\qquad R=m,\qquad U=m^2,\qquad U_1=m^3.
\]

No integer source \(N\) is constructed. In particular, these samples do
not certify \(A=\log N\), global CA optimality at the own price,
GA1, regularity, or the comparison with all larger integers.

With \(\rho_z(d)=\mu(d)\log(z/d)\mathbf1_{d\le z}\), the weights are

\[
\lambda=\rho_R/\log R,\qquad
\theta'=(\rho_{U_1}-\rho_U)/\log(U_1/U),\qquad
h(k)=\sum_{[d,e]=k}\lambda(d)\theta'(e).
\]

The measured quantities are the cumulative divisor response
\(H(y)=\sum_{n\le y}(\sum_{d\mid n}\lambda(d))
(\sum_{e\mid n}\theta'(e))\), its mean coefficient
\(S=\sum_k h(k)/k\), and the floor correction \(H(y)-yS\).
The program retains the floor factors across every cutoff.

The finite scans give the following numerical results. Each minimum
refers to all integer \(1\le y\le m^4\).

| \(m\) | Mean \(S\) | Minimum \(H(y)\) | Position \(y\) |
|---:|---:|---:|---:|
| 2 | 0.0727044767703 | 0.736965594166 | 6 |
| 3 | 0.0404523208756 | 0.738150879213 | 15 |
| 4 | 0.0333657291172 | 0.639694337482 | 39 |
| 6 | 0.0192798139577 | 0.288748998431 | 95 |
| 8 | 0.0132354556185 | 0.283280065293 | 222 |
| 12 | 0.00799243811993 | 0.639044295574 | 221 |
| 16 | 0.00561368643133 | 0.830482191778 | 346 |

The scans contain 92,017 integer points. No response below
\(-10^{-10}\) was observed. Selected values, including each grid
minimum and the larger arguments \(m^6,m^8\), are also evaluated by a
direct pair sum at 80 and 120 decimal digits. Their largest difference
is below \(7.3\cdot10^{-71}\). At the selected grid points, the direct
divisor-response computation and grouped LCM-floor computation differ
by at most \(1.14\cdot10^{-13}\).

These calculations supply neither a sign theorem for all real clocks
nor a bound on the actual-source weighted negative part. Decimal
agreement is not an interval enclosure or Lean certificate. Positive
mean values also do not determine the sign of the floor response at
unmeasured arguments. No asymptotic rate is fitted to these samples.

The published mean supplier is Carneiro–Chirre–Helfgott–Mejía-Cordero,
[arXiv:2005.03162v6](https://arxiv.org/abs/2005.03162v6), Theorem 1.1,
equations (1.7)–(1.8). Its first two terms cancel in the polarization
formula for this mixed mean. That supplier remains a classical result
to reuse; its two-sided subexponential remainder does not supply the
required same-source one-sided fixed-power bound.

Run from the experiment repository root:

```sh
python3 docs/reports/fib-robin-boundary/mixed_sieve_response.py
```

The JSON file contains the parameter values, support sizes, numerical
agreement checks, selected floor corrections and explicit scope
limits. The program requires only Python's standard library. It is an
original implementation of the finite definitions; no third-party
source code is copied.
