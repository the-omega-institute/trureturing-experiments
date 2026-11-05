# Cauchy conditional operator: exact checks and numerical diagnostics

The mathematical input is Olga V. Holtz, *The Marvelous Meixner Polynomials:
Second Kind, First Class*, arXiv:2611.01422v1, equations (5.1), (5.5), (14.1),
the conditional-moment formula (14.2), and Proposition 14.7. The source asks
for the essential spectrum in Section 16 without specifying a convention.
The proposed answer here uses failure of Fredholmness of Q − λI.

Run with Python 3, SymPy, NumPy and mpmath:

```sh
python3 check.py > check.out
python3 numerical.py > numerical.out
```

`check.py` uses exact arithmetic. It checks the conditional mean and variance,
the Cauchy interval-mass ratio 2/3, and the edge-gamma moments. For η = n/2,
1 ≤ n ≤ 8, it evaluates the Mellin integral for arbitrary symbolic real ξ and
obtains Γ(η+1/2+iξ)/(Γ(η)η^(1/2+iξ)). It also checks 56 nonnegative Mellin
moments, the logarithmic Jacobian, elementary remainder coefficients, the inverse-symbol algebra, and the cancellation and coefficient in the
log-Gamma remainder. The analytic remainder theorem itself is supplied by
Brent, arXiv:1609.03682v2, Section 2, Theorem 1; the script checks its application algebra.
It also checks two rational inequalities used to separate the Fredholm and
Weyl conventions, without certifying an operator index. Its final line is `ALL_CHECKS_OK`.
These identities do not establish the infinite-dimensional essential spectrum.

`numerical.py` uses floating point arithmetic solely as diagnostics. Its
Nyström matrices discretize the actual source density after conjugating to
Lebesgue L². The full sections include a core of 48 midpoint nodes and two
logarithmic tails, with 96 or 160 nodes each and tail lengths 6 or 10.
A separate positive-tail calculation compares the actual kernel with its
Wiener–Hopf limit and acts on sine-window Fourier packets at ξ = 0, 1/2, 1.
It uses step 1/12, shifts 2 and 4, and window lengths 8, 16 and 32.

| η | Proposed essential norm | Packet residual at ξ=0, shift=4, length=8 | Length=32 |
|---|---:|---:|---:|
| 1/2 | 0.797884560803 | 0.13077823 | 0.017361963 |
| 1 | 0.886226925453 | 0.091438257 | 0.010840657 |
| 3 | 0.959368788700 | 0.043629353 | 0.0050103213 |

At length 8, the finite-window Frobenius differences between the actual and
limiting tail matrices decrease from 0.0011238935, 0.0010600674, 0.0018104538
at shift 2 to 0.000020404338, 0.000019547537, 0.000033360551 at shift 4.
These are finite-window readings, not estimates for an infinite-region norm.
The largest full-section norm is 1.000019637; the source's true contraction
bound is 1, so the finite matrices retain quadrature error. Nonnormal
finite-section eigenvalues and floating point packets certify no Fredholm
or essential-spectrum assertion.

The proposed proof still needs a Hilbert–Schmidt compactness bridge, the
compact tail reduction, and a two-sided compressed Fourier parametrix.
No Lean proof of either the proposed answer or its negation is included.
The programs are original consistency checks and contain no copied
third-party program source.
