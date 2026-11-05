# Unconditional orthant projection

Exact arithmetic checks for Conjecture 3, equation (8), in Fradelizi–Manui–Meyer–Ndiaye, *L_p-Rogers–Shephard type inequalities for L_p-zonoids and symmetric bodies*, [arXiv:2607.03582v1](https://arxiv.org/abs/2607.03582v1).

Run with Python 3 and SymPy (checked with SymPy 1.14.0):

    python3 docs/reports/fmmn-unconditional-orthant-projection/check.py

`output.txt` is the successful stdout. The program uses exact integer/rational polynomial arithmetic. It checks endpoints and every stationary point of each polynomial of degree at most two, so its interval certificates cover every real parameter in [0,4]. It also cross-checks an independent rational convex hull implementation on 185 rational parameters. These are computational certificates, without Lean kernel verification.

For 0<t≤4, define B(g) as the box with coordinate radii g, and set

    g1 = (2,1,3,4), g2 = (4,4,0,2), g3 = (4−t,4,0,2+t),
    K_t = conv(B(g1) ∪ B(g2) ∪ B(g3)),
    Lx = (x1−x2+x4, x3+x4), E = row(L).

All boxes are invariant under coordinate sign changes; B(g1) has volume 384 and lies in K_t. Thus K_t is an unconditional full-dimensional convex body. The Gram determinant of L is 5, L=LP_E, and L restricted to E multiplies area by √5.

The orthant identity is K_t ∩ R⁴₊ = conv([0,g1] ∪ [0,g2] ∪ [0,g3]), where [0,g] is a coordinate box. To prove it, write a nonnegative x as Σλ_i v_i with v_i in B(g_i) and put s_j=Σλ_i|v_ij|. When s_j>0 replace v_ij by x_j|v_ij|/s_j; when s_j=0 use zero. Since 0≤x_j≤s_j, these replacement vectors belong to [0,g_i] and their weighted sum is x. The reverse inclusion follows from convexity.

The checker constructs all box corners, certifies their convex shadow polygons using edge halfplanes and vertex attainability, and computes areas by the shoelace formula. It independently reconstructs the eight-facet base polygon from rational halfspace intersections. The extra signed shadow adds two opposite triangles of area 3t/2 each, while its positive shadow stays inside the original pentagon. Consequently the positive area is 87/2 and the full area is 174+3t. The actual orthogonal projection ratio is 4+2t/29>4, with value 120/29 at t=2.

**Source fidelity:** Proposition 4's line upper bound is attained on the same family for the line spanned by (1,−1,0,1). The projected lengths are 20/√3 and 10/√3, respectively, giving ratio 2.

Cube padding gives ratio 2ⁿ(1+t/58) in ambient dimension n+2 for every n≥2. Taking k copies gives ratio 2²ᵏ(1+t/58)ᵏ in ambient dimension 4k. Projections and orthant intersections factor, and Euclidean product volume is multiplicative. For fixed t>0, Bernoulli's inequality gives (1+t/58)ᵏ≥1+kt/58, ruling out every fixed multiplier C in front of 2ⁿ. The program verifies the product identities for n=2,…,12 and k=1,…,12; the unbounded conclusion uses the general product and Bernoulli arguments.

The construction and certificates are research results. The program does not establish literature novelty or provide a Lean proof.
