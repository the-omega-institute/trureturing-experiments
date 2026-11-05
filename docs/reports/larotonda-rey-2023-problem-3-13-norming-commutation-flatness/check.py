#!/usr/bin/env python3
"""Exact checks for Larotonda--Rey Problem 3.13.

All evidence is in SymPy over exact integers/rationals/radicals.  The script
checks the trace identity used in the source and evaluates the proposed
Minkowski gauge on the n=2 counterexample.
"""

import sympy as sp


I = sp.I


def commutator(a, b):
    return a * b - b * a


def hermitian_trace_radius_squared(h):
    n = h.rows
    t = sp.trace(h) / n
    h0 = h - t * sp.eye(n)
    return sp.simplify(sp.trace(h0 * h0) / n)


def raw_gauge_from_t_and_r(t, r, a):
    """Gauge of t^2 + max(r-a,0)^2 <= 1 in the active branch r/lambda >= a."""
    # The chosen vectors have r >= a at the gauge boundary.  Solving
    # (t/lambda)^2 + (r/lambda-a)^2 = 1 gives this positive root.
    return sp.simplify(
        (sp.sqrt((1 - a * a) * t * t + r * r) - a * r) / (1 - a * a)
    )


def raw_gauge(h, a):
    n = h.rows
    t = sp.trace(h) / n
    r2 = hermitian_trace_radius_squared(h)
    r = sp.sqrt(r2)
    value = raw_gauge_from_t_and_r(t, r, a)
    # Confirm that the active-branch formula is applicable at lambda=value.
    assert sp.simplify(r / value - a) >= 0
    return value


def check_source_trace_identity():
    # This is the cyclic identity in Definition 2.1's trace-duality paragraph:
    # (V | [X,Y]) = (Y | [V,X]) for skew-Hermitian matrices.
    x = sp.Matrix([[0, 1], [-1, 0]])
    y = I * sp.Matrix([[1, 2], [2, -1]])
    v = I * sp.Matrix([[2, -1], [-1, 3]])
    lhs = sp.simplify(-sp.trace(v * commutator(x, y)))
    rhs = sp.simplify(-sp.trace(y * commutator(v, x)))
    assert lhs == rhs
    return lhs


def check_frobenius_equality_fidelity():
    # Theorem 3.11 / Corollary 3.12 specialize to the Frobenius norm: a
    # nonzero commutator moves a vector strictly off the equality case.
    k = sp.diag(1, -1)
    x_commuting = sp.diag(2, -3)
    x_noncommuting = sp.Matrix([[0, 1], [-1, 0]])
    v = I * k

    def frob_sq(z):
        return sp.simplify(sp.trace(z.conjugate().T * z))

    same = sp.simplify(frob_sq(v + commutator(x_commuting, v)) - frob_sq(v))
    moved = sp.simplify(frob_sq(v + commutator(x_noncommuting, v)) - frob_sq(v))
    assert same == 0
    assert moved > 0
    return same, moved


def check_counterexample():
    n = 2
    a = sp.Rational(1, 2)
    k = sp.diag(1, -1)
    x = sp.Matrix([[0, 1], [-1, 0]])
    v_h = sp.eye(n) + a * k
    v = I * v_h
    displacement = commutator(x, v)

    # K has simple spectrum, r(K)=1 and 0<a<1/||K||op exactly.
    assert hermitian_trace_radius_squared(k) == 1
    assert a > 0 and a < 1
    assert v_h[0, 0] != v_h[1, 1]

    # The unique support at V is central N=iI/n.  It commutes with every X_C.
    norming_n = I * sp.eye(n) / n
    assert commutator(norming_n, x) == sp.zeros(n)
    assert sp.simplify(-sp.trace(norming_n * v)) == 1

    f_v_raw = raw_gauge(v_h, a)
    assert f_v_raw == 1

    # More directly, v+s[X,V] = i times this Hermitian part; use a symbolic
    # parameter for the full straight commutator line.
    s = sp.symbols("s", real=True)
    line_h = v_h + (-I * s * displacement)
    r2_line = sp.factor(hermitian_trace_radius_squared(line_h))
    # Here Tr([X,K]^2)/n = 4, so the source's q-factor is 4.
    assert sp.simplify(r2_line - a * a * (1 + 4 * s * s)) == 0

    # At every nonzero s the boundary expression is strictly above one.
    # Its exact value is 1 + (a*(sqrt(1+s^2)-1))^2.
    excess = sp.simplify(
        (sp.sqrt(r2_line) - a) ** 2
    )
    assert excess.subs(s, 1) > 0
    assert raw_gauge(line_h.subs(s, 1), a) > f_v_raw

    # Positive rescaling enforces the source's rank-one normalization ||iP||=1.
    p = sp.diag(1, 0)
    rank_one_raw = raw_gauge(p, a)
    scale = sp.simplify(1 / rank_one_raw)
    assert rank_one_raw == (sp.sqrt(7) - 1) / 3
    assert scale > 0
    assert sp.simplify(scale * raw_gauge(v_h, a)) == scale
    assert sp.simplify(scale * raw_gauge(line_h.subs(s, 1), a)) > scale

    return {
        "raw_norm_V": f_v_raw,
        "raw_norm_rank_one": rank_one_raw,
        "normalization_scale": scale,
        "radius_squared_on_line": r2_line,
        "strict_excess_at_s1": excess.subs(s, 1),
        "raw_norm_at_s1": raw_gauge(line_h.subs(s, 1), a),
        "commutator": displacement,
    }


def main():
    trace_value = check_source_trace_identity()
    frob_same, frob_moved = check_frobenius_equality_fidelity()
    witness = check_counterexample()

    print("SOURCE_TRACE_IDENTITY", trace_value)
    print("FROBENIUS_EQUALITY_DIFFERENCES", frob_same, frob_moved)
    for key, value in witness.items():
        print(key.upper(), value)
    print("ALL_CHECKS_OK")


if __name__ == "__main__":
    main()
