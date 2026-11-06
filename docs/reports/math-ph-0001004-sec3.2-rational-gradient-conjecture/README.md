# Rational-gradient first integrals: exact candidate checks

Run `python3 check.py` with Python 3 and SymPy. `check.out` contains its output.
The program uses exact symbolic algebra and integers. It neither proves the
conjecture nor its negation in Lean.

Source: L. G. S. Duarte, L. A. da Mota and J. E. F. Skea, arXiv:math-ph/0001004v1,
Section 3.2, equations (25)–(26). The conjecture asks whether two independent
nonconstant elementary first integrals of a rational scalar second-order ODE
force a nonconstant first integral r₀ + Σ cᵢ log rᵢ with rational rᵢ.
The associated journal article has an additional author, S. E. S. Duarte,
DOI 10.1088/0305-4470/34/14/308. These checks concern the arXiv statement.

## Exact results and limits

Fidelity: the source's harmonic-oscillator example has y″ = −y, S = y/y′,
R = y′ and the conserved energy C1 = y² + y′². The program verifies its three
compatibility equations and energy conservation. A phase invariant is an
additional check, not an additional result attributed to the source.

For symbolic p,q, the program verifies the counterexample's coordinate
transformations, Riccati equations, jet-coordinate Jacobian and inverse,
and the complete displayed rational right-hand side. It verifies both
logarithmic first integrals and the nonzero coefficient 1/(vw) of their
exterior product. It verifies the general-solution Riccati identity, the
closed-horizontal differential identity and inversion's action on dlog.
The irreducible polynomial T⁴ − 10T² + 1 supplies an exact degree-four
number-field check. The finite monomial-weight checks are diagnostics;
they do not establish universal nonexistence.

The universal obstruction is the mathematical argument below. SymPy checks
its differential identities, not a quantification over all rational functions.
No kernel verification or certified literature novelty is asserted.

## Counterexample and universal obstruction

Let m ≥ 1 be an integer, α = √2, β = m√3, p = 2 and q = 3m². The numbers
1,α,β are linearly independent over Q: the biquadratic field
Q(√2,√3) has basis 1,√2,√3,√6, and scaling √3 by nonzero m preserves
independence. The degree-four minimal polynomial of √2 + √3 is checked
exactly by the program.

Put F = C(u,v,w), E = u∂u + αv∂v + βw∂w, and let σ invert u,v,w.
Set s = (u−1)/(u+1), r = (v−1)/(v+1), t = (w−1)/(w+1),
x = s², A = αr/s and B = βt/s. Then F = C(s,A,B), σ(s) = −s,
and σ fixes A,B. Consequently F^σ = K = C(x,A,B), a rational field,
and F/K is quadratic. The inverse formulas are
u = (1+s)/(1−s), v = (α+sA)/(α−sA), w = (β+sB)/(β−sB).

Since E(x) = s(1−x) and σ reverses both E and s(1−x),
D = E/[s(1−x)] descends to K and D(x) = 1. With

- a = −1/[2(1−x)], b = −1/(2x), c_h = h/[2x(1−x)],
- D(A) = aA² + bA + c_p and D(B) = aB² + bB + c_q,
- y = A−B, z = D(y), H = A+B = (z−by−c_p+c_q)/(ay),

we have K = C(x,y,z). The Jacobian of (A,B) ↦ (y,z) is 2ay, nonzero
as a rational function. Differentiating z gives

Phi = 3z²/(2y) + [(p−q)/(xy)−1]z/(x−1) + y³/[8(x−1)²]
+ [3−2(p+q+1)x−x²]y/[8x²(x−1)²] + (p−q)/[2x²(x−1)]
+ (p−q)²/[8x²y(x−1)²].

This expression has rational coefficients for every such m. At m = 1,
x = 1/4, A = 2, B = 1 gives (y,z,Phi) = (1,−20/3,328/9),
with the displayed chart and logarithmic differential denominators nonzero.

**Rational constants.** Suppose f = P/Q ∈ F, P,Q coprime and nonzero,
and E(f) = 0. Then E(P)Q = PE(Q), so unique factorization gives
P | E(P) and Q | E(Q). The operator preserves or lowers total degree;
therefore E(P) = λP and E(Q) = λQ for λ ∈ C, allowing λ = 0.
The weights i + αj + βk of distinct polynomial monomials are distinct.
Every eigenpolynomial is therefore one monomial. Equal weights for P,Q
force their exponent triples to agree; f is constant. The zero function
is already constant. Thus the rational constant field of E is C.

**Closed horizontal forms.** Write η = Σ hᵢ dUᵢ, with U = (u,v,w),
λ = (1,α,β), dη = 0 and η(E) = Σ λᵢUᵢhᵢ = 0. Differentiating the last
identity with respect to Uⱼ and using ∂ⱼhᵢ = ∂ᵢhⱼ gives
E(hⱼ) + λⱼhⱼ = 0. Hence E(Uⱼhⱼ) = 0. The preceding result gives
Uⱼhⱼ = cⱼ ∈ C, and η = Σ cⱼ dlog Uⱼ, with Σ λⱼcⱼ = 0.

**Descent obstruction.** If a nonzero closed D-horizontal rational form
existed over K, its pullback to F would be one of the forms just found.
Every such form is negated by σ, since σ*(dlog Uᵢ) = −dlog Uᵢ.
A form pulled back from K is σ-invariant. Characteristic zero then forces
it to vanish. Pullback is injective for the separable algebraic extension
F/K (equivalently use the bases ds,dA,dB and dx,dA,dB, with dx = 2s ds).
Thus no nonzero closed rational D-horizontal form exists over K.

**Elementary integrals.** In a local meromorphic function field, choose
branches of s = √x and logarithms on a connected regular chart. The field
obtained from K by adjoining s and log u, log v, log w is elementary for
the commuting coordinate partial derivations. Its common constants are C:
a meromorphic germ killed by all coordinate partials is a constant germ.
The two functions I₁ = log v − α log u and I₂ = log w − β log u satisfy
D(Iᵢ) = 0 and have independent differentials; the dv∧dw coefficient is
1/(vw). They therefore meet the conjecture's hypotheses. Any nonconstant
J = r₀ + Σ cᵢ log rᵢ with rational rᵢ would give a nonzero closed rational
form dJ annihilating D, contradicting the descent obstruction.

The full elementary solution is, for s = √x and U = (1+s)/(1−s),
y = α/s (C₁ U^α−1)/(C₁ U^α+1) − β/s (C₂ U^β−1)/(C₂ U^β+1),
where U^γ = exp(γ log U). Its two parameters vary A and B independently.

These arguments concern the old rational-gradient conjecture. They do not
refute rational-S algorithms or methods allowing algebraic integrating factors.
The general solution and all field formulas are local, away from their poles.
