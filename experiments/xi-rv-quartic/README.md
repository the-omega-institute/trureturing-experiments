# Original-theta quartic RV numerator

`certify.py` evaluates the actual profile in Silva,
[arXiv:2609.25564v1](https://arxiv.org/html/2609.25564v1), equations
(9), (10) and (16), including both the infinite theta-series remainder
and the unbounded integration tail. It reuses the existing theta-jet
and derivative-envelope suppliers under `docs/reports/theta-mixed-matrix`.
No finite theta truncation is substituted for the source.

Run from this repository root:

```sh
uv run --no-project --python python3.12 --with python-flint==0.9.0 python experiments/xi-rv-quartic/certify.py --output experiments/xi-rv-quartic/certificate.json
```

With `a=P(1/4)`, `b=P(1/2)` and `P(0)=1/2`, the reciprocal numerator is

$$
U_4(q)=\tfrac12(1+q^4)-4a(q+q^3)+6bq^2.
$$

For `q` on the unit circle, `t=q+q^{-1}` is real and

$$
q^{-2}U_4(q)=\tfrac12(t^2-2)-4at+6b.
$$

The output supplies strict rational enclosures with

$$
16a^2-12b+2<0.
$$

Thus this real quadratic has no real zero, and `U4` has no unit-circle
zero. Its constant term is nonzero, so all four numerator zeros lie
off the unit circle. The certificate also checks `U4(1)>0`.

The universal claim that every even-degree numerator has all zeros
on the unit circle fails at E=4. This does not exclude a sufficiently
large-degree subsequence, prove that the transformed polynomial has an off-line
zero, or prove or refute RH. Arb enclosures are numerical certificates;
no Lean verification is supplied.
