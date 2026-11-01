# Yang–Morales bounded shadowing converse

Source: Y. Yang and C. A. Morales, *Topological shadowing for linear operators*,
[arXiv:2608.12862v1](https://arxiv.org/abs/2608.12862v1), Definition 3 and the
question immediately following it. The candidate answer is negative: for every
real pair $0<a<1<b$, the complex operator $T(u,v)=(au,bv)$ on
$\mathbb C^2$ with the maximum norm has topological bounded shadowing and
fails topological shadowing.

Run with Python 3 and its standard library:

```sh
python3 docs/reports/yang-morales-bounded-shadowing-converse/check.py
```

The stored output is `check.out`. All evidence arithmetic uses integers and
`fractions.Fraction`.

- Fidelity to the proof of source Theorem A: reproduce its finite
  variation-of-constants formula
  $x_n=L^nx_0+\sum_{j=0}^{n-1}L^jr_{n-j}$ and its inverse-normalized
  version $L^{-n}x_n=x_0+\sum_{i=1}^nL^{-i}r_i$ for $L(z)=2z$ on
  $\mathbb C$, every five-letter error word over $\{-1/7,0,1/7\}$,
  and times 1 through 5: 1,215 exact checks of both displayed identities.
- Fidelity to source Lemma 5: four complex rotations (multipliers $1,-1,i,-i$)
  admit the radius-3 ramp obstruction with strict errors for the continuous
  tolerance $\delta(z)=1/(1+|z|^2)$. The ramp has 31 steps. Its endpoint
  tracing balls of radius 1 are incompatible. This checks concrete instances
  of the lemma's mechanism, not its quantifier over every positive tolerance.
- Bounded recurrence: 132 nondegenerate rational parameter pairs and all
  5,148 coordinate error words of periods 1 through 3 over
  $\{-\delta/2,0,\delta/2\}$ satisfy the exact cyclic Green formulas.
  The checks include both recurrence directions, four endpoint lengths,
  and the zero-tracer bound $\|x_n\|_\infty\le 1/2$ for $\epsilon(0)=1$.
- Full-shadowing obstruction: 396 rational splice certificates use three
  nonconstant control functions, evaluated at $Tx_n$, and tracing tolerance
  $\epsilon(z)=1/(1+\|z\|_\infty^k)$, evaluated on the forced tracer's orbit.
  An exact integer failure time is found in each case; the largest is 391.

The mathematical argument additionally uses convergence of the endpoint terms
for arbitrary bounded two-sided sequences, continuity of arbitrary tolerances,
and divergence of geometric powers. These universal steps are described in the
preregistration; the finite checks do not prove them. No Lean proof or complete
finite-dimensional classification is supplied by this experiment.
