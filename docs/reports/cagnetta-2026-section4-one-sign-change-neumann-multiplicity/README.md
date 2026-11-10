# One sign change and three positive Neumann solutions

`check.py` uses only Python integers and `fractions.Fraction`. Run it with Python
3.12 or newer. Its captured standard output is `check.out`; success ends with
`ALL_CHECKS_OK`.

The source is Alberto Cagnetta, *A uniqueness result for two-step indefinite
p-Laplacian equations*, arXiv:2610.06068v1, Section 4. The fidelity checks reproduce
its weight integral −35877/320 and its two positive-hump maxima 13/800 and 301/80.
These checks do not certify the source's three numerical shooting roots.

The independent construction takes γ=2, b=10000, r=1/b, θ=7/8 and
C=b/θ²=640000/49. Put s=1−z and D=z−(1−z)(z²+r). Composing a kick of mass 1,
a drift of length 1, a kick of mass b, a drift of length r and a kick of mass −C
produces terminal velocity

    F₀(z) = s² [C D² − (1+b z²)].

The program checks this composition directly and certifies alternating strict
signs at z=7/1600, 1/8, 1/2, 7/8. The polynomial identity

    D − (3z/4−r) = z(z−1/2)² + rz

proves D≥509/160000 throughout that interval. Consequently all ideal trajectories
have u≥δ=509/1280000, not just the four sign-testing trajectories.

## Finite coefficients

Set ε=2^(−650000), T=1+r+ε. On the consecutive intervals of lengths
ε, 1−ε, ε, r−ε, ε, take the coefficient values

    1/ε+ε, ε, b/ε+ε, −ε, −C/ε−ε.

All lengths are positive, the coefficient changes sign exactly once, and its
integral is 1+b−C+ε(1+ε−r)<0. These are finite, bounded coefficients.

## Analytic bridge certified by the arithmetic

Parametrize each stage by τ∈[0,1], with state (u,v). The actual fields are

    (εv, −(1+ε²)u²)
    ((1−ε)v, −(1−ε)εu²)
    (εv, −(b+ε²)u²)
    ((r−ε)v, +(r−ε)εu²)
    (εv, +(C+ε²)u²).

The corresponding ideal fields are (0,−u²), (v,0), (0,−bu²),
(rv,0), (0,+Cu²). On δ/2≤u≤2 and |v|≤30000 their sup-norm Lipschitz
constants are at most L=60000; actual-minus-ideal fields have norm at most 30000ε.
The program checks the inequalities 4(C+1)<L, 1+b+C<29999 and all remaining
rational margins. Since e<4, the five-stage Gronwall recurrence

    E₀=0; Eⱼ₊₁=(Eⱼ+30000ε) 2^120000

bounds every stage's state error, and E₅≤150000·2^(−50000). This is less than
δ/2, 1 and half each of the four nonzero ideal terminal velocities.

Local existence for polynomial vector fields, followed by this estimate up to
the first exit, prevents exit from the rectangle and gives continuation through
all five stages. Continuous dependence makes terminal velocity continuous in z.
The four signs therefore persist for the actual finite coefficients. The
intermediate value theorem gives a zero in each of the three intervening
intervals. Their distinct initial values 1−z give three distinct solutions.
Their u is positive, u is C¹ and v=u′ is absolutely continuous. Rescaling
U(t)=u(Tt), A(t)=T²a(Tt) gives the normalization T=1.

The program checks exact algebra and all numerical hypotheses of this analytic
argument. Picard existence, continuation, Gronwall and the intermediate value
theorem are mathematical dependencies; they are not machine-checked Lean
conclusions of this report. No floating-point integration is used as evidence.
