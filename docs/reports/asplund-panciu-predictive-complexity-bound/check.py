#!/usr/bin/env python3
"""Exact source-fidelity and proposed-refutation checks for arXiv:2309.15200v4.

Requires Python 3 and SymPy. All evidence is integer/rational/symbolic algebra.
No floating point evaluation or Lean proof of the conjecture is used.
"""

import sys
from fractions import Fraction

import sympy as sp


def require(condition, label):
    if not condition:
        raise AssertionError(label)


def same(left, right, label):
    difference = left - right
    entries = list(difference) if isinstance(difference, sp.MatrixBase) else [difference]
    require(all(sp.simplify(x) == 0 for x in entries), label)


def norm2(z):
    return sp.simplify(z * sp.conjugate(z))


def source_merge(rows):
    """Equations (18)/(23): merge columns 1 and 2 into beta, keep column 3."""
    result = []
    for row in rows.tolist():
        phase_sum = sp.simplify(row[0] + row[1])
        mass = sp.simplify(norm2(row[0]) + norm2(row[1]))
        require(phase_sum != 0, "phase denominator is nonzero")
        coefficient = sp.sqrt(mass) * phase_sum / sp.sqrt(norm2(phase_sum))
        result.append([coefficient, row[2]])
    return sp.Matrix(result)


def marginal(rows):
    """Independent row-Gram realization of the literal partial trace."""
    return sp.simplify(rows * rows.adjoint())


def direct_partial_trace(vector, a, b):
    joint = vector * vector.adjoint()
    return sp.Matrix(a, a, lambda i, k: sum(joint[i * b + j, k * b + j]
                                          for j in range(b)))


def fidelity_check():
    # Source Section 2.4, immediately after (26): Bell state 11 + 22.
    rows = sp.Matrix([[1 / sp.sqrt(2), 0, 0], [0, 1 / sp.sqrt(2), 0]])
    merged = source_merge(rows)
    original_rho, merged_rho = marginal(rows), marginal(merged)
    same(original_rho, sp.eye(2) / 2, "source Bell marginal")
    same(merged_rho, sp.ones(2) / 2, "source merged Bell marginal")
    require(original_rho.eigenvals() == {sp.Rational(1, 2): 2}, "Bell spectrum")
    require(merged_rho.eigenvals() == {sp.Integer(0): 1, sp.Integer(1): 1},
            "source predictive Bell spectrum")
    print("FIDELITY: Section 2.4 Bell example: spectra {1/2,1/2} -> {1,0}; S_A=log(2), C_A=0")


def dynamics_check():
    # z=exp(-i*tau). Algebra is in Laurent expressions, with conjugate(z)=1/z.
    z = sp.Symbol("z", complex=True, zero=False)
    x0, x1, b0, b1, b2, c0, c1, c2 = sp.symbols("x0 x1 b0 b1 b2 c0 c1 c2", complex=True)
    rho = sp.Matrix([x0, x1]) * sp.Matrix([x0, x1]).adjoint()
    propagator = sp.diag(1, 1, z, 1, 1, 1 / z)
    adj = propagator.adjoint().subs(sp.conjugate(z), 1 / z)
    same(propagator * adj, sp.eye(6), "propagator unitarity")
    ket = sp.Matrix([x0 * b0, x0 * b1, x0 * b2, x1 * b0, x1 * b1, x1 * b2])
    evolved = propagator * ket
    reduced = direct_partial_trace(evolved, 2, 3).subs(sp.conjugate(z), 1 / z)
    d = sp.diag(z, 1 / z)
    dadj = sp.diag(1 / z, z)
    expected = (norm2(b0) + norm2(b1)) * rho + norm2(b2) * d * rho * dadj
    same(reduced, expected, "exact channel: depends only on p=|b2|^2 for unit b")
    v = sp.Matrix([x0 * b0, x0 * b1, 0, x1 * b0, x1 * b1, 0])
    w = sp.Matrix([x0 * c0, x0 * c1, 0, x1 * c0, x1 * c1, 0])
    cross = propagator * v * w.adjoint() * adj
    traced = sp.Matrix(2, 2, lambda i, k: sum(cross[3 * i + j, 3 * k + j]
                                             for j in range(3)))
    overlap = b0 * sp.conjugate(c0) + b1 * sp.conjugate(c1)
    same(traced, overlap * rho, "Equation (8) for orthogonal vectors in W")
    general_w = sp.Matrix([x0 * c0, x0 * c1, x0 * c2,
                           x1 * c0, x1 * c1, x1 * c2])
    general_cross = propagator * ket * general_w.adjoint() * adj
    general_traced = sp.Matrix(2, 2, lambda i, k:
                              sum(general_cross[3 * i + j, 3 * k + j]
                                  for j in range(3)))
    same(general_traced, overlap * rho + b2 * sp.conjugate(c2) * d * rho * dadj,
         "general cross trace: no additional coherent equivalent orthogonal pair")
    p, q = sp.symbols("p q", real=True)
    plus = sp.ones(2) / 2
    channel = lambda weight: (1 - weight) * plus + weight * d * plus * dadj
    same((channel(p) - channel(q))[0, 1], (p - q) * (z**2 - 1) / 2,
         "distinct p are distinguished whenever z^2 != 1")
    # Interacting dynamics: evolve a product with A=(1,1), B=(1,0,1).
    # The 2x2 minor in columns 1 and 3 is 1/z-z, nonzero at z=-i.
    entangling_rows = sp.Matrix([[1, 0, z], [1, 0, 1 / z]])
    require(entangling_rows[:, [0, 2]].det().subs(z, -sp.I) != 0,
            "Hamiltonian is genuinely interacting")
    print("DYNAMICS: H=Z tensor |3><3|; exact p-channel and Equation (8) verified symbolically")
    print("DYNAMICS: distinct p distinguished by (p-q)(z^2-1)/2; a product input becomes entangled")


def family_check():
    u, v, t, gap = sp.symbols("u v t gap", positive=True)
    s, c = u**2 + v**2, 2 * u * v
    denominator = 2 * (s + t**2)
    rows = sp.Matrix([[u, v, t], [v, u, -t]]) / sp.sqrt(denominator)
    merged = source_merge(rows)
    rho, rho_prime = marginal(rows), marginal(merged)
    expected = sp.Matrix([[sp.Rational(1, 2), (c - t**2) / denominator],
                          [(c - t**2) / denominator, sp.Rational(1, 2)]])
    expected_prime = sp.Matrix([[sp.Rational(1, 2), (s - t**2) / denominator],
                                [(s - t**2) / denominator, sp.Rational(1, 2)]])
    same(rho, expected, "whole original family")
    same(rho_prime, expected_prime, "literal source-merge family")
    same(direct_partial_trace(sp.Matrix(list(rows)), 2, 3), rho, "independent partial trace")
    same(direct_partial_trace(sp.Matrix(list(merged)), 2, 2), rho_prime,
         "independent predictive partial trace")
    same(sp.trace(rho), 1, "original normalization")
    same(sp.trace(rho_prime), 1, "merged normalization")
    factor = (u - v)**2 * (2 * t**2 - (u + v)**2) / (4 * (s + t**2)**2)
    same(rho_prime.det() - rho.det(), factor, "determinant gain factorization")
    same(((s - t**2)**2 - (c - t**2)**2),
         (u - v)**2 * ((u + v)**2 - 2 * t**2), "squared-coherence factorization")
    # Under t^2=(u+v)^2/2+gap and u!=v, gain/(u-v)^2 is strictly positive.
    positive_gain = sp.factor((factor / (u - v)**2).subs(t**2, (u + v)**2 / 2 + gap))
    require(positive_gain.is_positive is True, "open-family strict sign")
    # Full rank follows without eigenvalue approximations.
    same(4 * (s + t**2)**2 * rho.det(), (u + v)**2 * ((u - v)**2 + 2 * t**2),
         "original determinant positive for u!=v")
    same((s + t**2)**2 * rho_prime.det(), s * t**2, "merged determinant positive")
    p = sp.Symbol("p", positive=True)
    h = -p * sp.log(p) - (1 - p) * sp.log(1 - p)
    same(sp.diff(h, p), sp.log(1 - p) - sp.log(p), "binary entropy derivative")
    # For p in (1/2,1), 0<1-p<p, hence this derivative is strictly negative.
    # Thus exact strict reduction in |offdiagonal| implies a strict entropy increase.
    print("FAMILY: det(rho_prime)-det(rho)=(u-v)^2*(2*t^2-(u+v)^2)/(4*(u^2+v^2+t^2)^2)")
    print("FAMILY: strict gain for u,v,t>0, u!=v, t^2>(u+v)^2/2; both marginals full rank")
    print("ENTROPY: h'(p)=log(1-p)-log(p)<0 on (1/2,1); more balanced spectrum has larger entropy")


def exact_witness_and_grid():
    rows = sp.Matrix([[3, 1, 3], [1, 3, -3]]) / sp.sqrt(38)
    rho, rho_prime = marginal(rows), marginal(source_merge(rows))
    require(sorted(rho.eigenvals()) == [sp.Rational(8, 19), sp.Rational(11, 19)],
            "exact original witness spectrum")
    require(sorted(rho_prime.eigenvals()) == [sp.Rational(9, 19), sp.Rational(10, 19)],
            "exact merged witness spectrum")
    # 19*(C-S)=log((11^11*8^8)/(10^10*9^9)); positivity is an integer comparison.
    numerator, denominator = 11**11 * 8**8, 10**10 * 9**9
    require(numerator > denominator, "strict entropy increase by integer powers")
    print("WITNESS: (u,v,t)=(3,1,3); spectra {8/19,11/19} -> {9/19,10/19}")
    print(f"EXACT_ENTROPY: exp(19*(C-S))={numerator}/{denominator}>1")
    count = 0
    for u in range(1, 9):
        for v in range(1, 9):
            for t in range(1, 13):
                if u == v or 2 * t * t <= (u + v)**2:
                    continue
                d = 2 * (u * u + v * v + t * t)
                r = Fraction(2 * u * v - t * t, d)
                r_prime = Fraction(u * u + v * v - t * t, d)
                require(abs(r_prime) < abs(r) < Fraction(1, 2), "exact rational grid")
                count += 1
    print(f"RATIONAL_GRID: {count} admissible integer triples checked (u,v=1..8; t=1..12)")


def main():
    print(f"Python {sys.version.split()[0]}; SymPy {sp.__version__}; exact arithmetic only")
    fidelity_check()
    dynamics_check()
    family_check()
    exact_witness_and_grid()
    print("ALL_CHECKS_OK")


if __name__ == "__main__":
    main()
