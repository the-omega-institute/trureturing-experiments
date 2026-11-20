"""Exact checks for the LorEA₂ convexity counterexample in arXiv:2506.14480v1."""

from itertools import product

import sympy as sp

Q = sp.Rational
I2 = sp.eye(2)
HADAMARDS = tuple(
    sp.Matrix(2, 2, [(-1 if i == j else 1) for i in range(4)])
    for j in range(4)
)
PHI1 = sp.diag(1, Q(1, 5))
PHI2 = sp.Matrix([[1, Q(3, 10)], [0, Q(1, 2)]])
R = sp.Matrix([[12, 5], [5, -12]]) / 13
H = HADAMARDS[-1]


def alpha_fourth(matrix):
    plus, minus = [matrix * sp.Matrix([1, sign]) for sign in (1, -1)]
    products = [column.dot(column) ** 2 for column in (plus, minus)]
    # For a real 2 × 2 matrix B, (σ₁ + σ₂)² = tr(BᵀB) + 2 |det B|.
    blocks = [matrix * hadamard * matrix.T for hadamard in HADAMARDS]
    entangled = [(b.T * b).trace() + 2 * abs(b.det()) for b in blocks]
    return max(products), max(entangled)


def check_fidelity():
    a1 = max(alpha_fourth(PHI1))
    a2 = max(alpha_fourth(PHI2))
    a12 = max(alpha_fourth(PHI1 + PHI2))
    # Separate the fourth roots by strictly ordered rational bounds.
    upper1, upper2, lower12 = Q(1056, 1000), Q(1394, 1000), Q(246, 100)
    assert a1 < upper1**4 and a2 < upper2**4
    assert lower12**4 < a12 and upper1 + upper2 < lower12
    print(f"Appendix B fidelity: alpha^4 = {a1}, {a2}, and {a12}")
    print(f"Appendix B strict inequality: alpha(phi1) + alpha(phi2) < "
          f"{upper1 + upper2} < {lower12} < alpha(phi1 + phi2)")


def check_base_case():
    f1, f2 = Q(50, 53) * PHI1, Q(5, 7) * PHI2
    assert alpha_fourth(PHI1) == (Q(676, 625), Q(776, 625))
    assert alpha_fourth(PHI2) == (Q(9409, 2500), Q(8969, 2500))
    for i, (matrix, alpha4, norm2, det) in enumerate([
        (f1, Q(7760000, 7890481), Q(2600, 2809), Q(500, 2809)),
        (f2, Q(9409, 9604), Q(97, 98), Q(25, 98)),
    ], 1):
        positivity4, entangled4 = alpha_fourth(matrix)
        assert max(positivity4, entangled4) == alpha4 < 1
        assert positivity4 == norm2**2 and 0 < norm2 < 1
        assert sp.diag(1, matrix.T).det() == det > 0
        print(f"Q{i}: alpha^4 = {alpha4}, positive-map norm^2 = {norm2}, det = {det}")

    # Identify a tensor with its coefficient matrix in the displayed basis.
    # The Lorentz generator w = (id ⊗ (1 ⊕ -R)) Î₂ is symmetric.
    assert R.T * R == I2 and R.T == R
    witness = sp.diag(1, -R)
    separator = sp.diag(1, H)
    diamond_rays = [sp.Matrix([1, sign * (i == 0), sign * (i == 1)])
                    for i, sign in product(range(2), (1, -1))]
    ray_values = {(x.T * separator * y)[0] for x, y in product(diamond_rays, repeat=2)}
    assert ray_values == {0, 2}
    print(f"R orthogonal and D₂ product-ray separator values are {sorted(ray_values)}")

    t, a, b = sp.symbols("t a b", real=True)
    phi = a * PHI1 + b * PHI2
    trace = sp.expand((R * phi * H * phi.T).trace())
    assert sp.expand(trace - (Q(362, 25)*a*a + 40*a*b + Q(1231, 50)*b*b)/13) == 0
    mixture = sp.diag(1, (t * f1 + (1 - t) * f2).T)
    paired = sp.expand((separator.T * mixture * witness * mixture.T).trace())
    assert sp.expand(paired - (1 - trace.subs({a: Q(50, 53)*t, b: Q(5, 7)*(1-t)}))) == 0
    assert sp.expand(paired - (414521*t*t - 504242*t + 120787)/3578666) == 0
    endpoints = (Q(1, 3), Q(5, 6))
    values = [paired.subs(t, u) for u in (*endpoints, Q(1, 2))]
    assert values == [Q(-5561, 16103997), Q(-415903, 128831976), Q(-110815, 14314664)]
    assert sp.Poly(paired, t).LC() > 0 and all(v < 0 for v in values)
    # Exact interpolation error proves the sign over the entire interval.
    left, right = endpoints
    chord = ((right - t)*values[0] + (t - left)*values[1]) / (right - left)
    assert sp.expand(chord - paired -
                     sp.Poly(paired, t).LC()*(t-left)*(right-t)) == 0
    print(f"separating polynomial = {sp.Poly(paired, t).as_expr()}")
    print(f"separator at 1/3, 5/6 and 1/2 = {values}")
    return f1.T, f2.T, paired


def check_retractions(spatial1, spatial2, paired):
    for n, m in product(range(2, 6), repeat=2):
        # pi drops coordinates, and its transpose inserts zeros.
        pi_input, pi_output = sp.eye(n + 1)[:3, :], sp.eye(m + 1)[:3, :]
        assert pi_input * pi_input.T == sp.eye(3)
        assert pi_output * pi_output.T == sp.eye(3)
        t = next(iter(paired.free_symbols))
        for spatial in (spatial1, spatial2, t * spatial1 + (1 - t) * spatial2):
            base = sp.diag(1, spatial)
            extended = pi_output.T * base * pi_input
            assert pi_output * extended * pi_input.T == base
            assert extended[0, :] == sp.eye(n + 1)[0, :]
        witness = pi_input.T * sp.diag(1, -R) * pi_input
        separator = pi_output.T * sp.diag(1, H) * pi_output
        mixture = pi_output.T * sp.diag(1, t*spatial1 + (1-t)*spatial2) * pi_input
        assert sp.expand((separator.T * mixture * witness * mixture.T).trace() - paired) == 0
    print("coordinate retractions and separating polynomial checked for 2 ≤ n,m ≤ 5")


def main():
    check_fidelity()
    check_retractions(*check_base_case())
    print("ALL_CHECKS_OK")


if __name__ == "__main__":
    main()
