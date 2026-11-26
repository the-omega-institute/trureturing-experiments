"""Finite certificates for the five-step overlap cones and their sharp free gap."""
from collections import Counter
from itertools import product
from math import comb
import sympy as sp

Q = sp.Rational
root5 = sp.sqrt(5)
lower_mass = (25 - 5*root5)/4
upper_mass = Q(7, 2)
high_mass = (25 + 5*root5)/4


def cone_table():
    # Every root from the minority-branch equation has order dividing 30.
    # At order 30, conjugation followed by negation sends a to 15-a.
    def sine_class(a):
        return min(a, (15-a) % 30)

    classes = [sine_class(a) for a in range(30)]
    cosines = {0: 1, 6: (root5-1)/4, 12: -(1+root5)/4,
               5: Q(1, 2), 15: -1}
    steps = sp.eye(4).col_join(-sp.ones(1, 4))
    gram = sp.eye(4) - sp.ones(4)/5
    assert steps*gram*steps.T == sp.eye(5) - sp.ones(5)/5
    assert steps.T*steps == gram.inv()
    table = Counter()
    light = set()
    for p in product(range(30), repeat=4):
        phases = (*p, -sum(p) % 30)
        if any(classes[a] != classes[p[0]] for a in phases[1:]):
            continue
        major, frequency = Counter(phases).most_common(1)[0]
        n = 5-frequency
        cos_major = cosines[min(major, 30-major)]
        signs = [1 if a == major else -1 for a in phases]
        jacobian = steps.T*sp.diag(*signs)*steps
        assert jacobian.det() == (-1)**n*(5-2*n)
        mass = sp.expand(5 - (5-2*n)*cos_major)
        orientation = (-1)**n
        table[mass, orientation] += 1
        if mass < Q(139, 40):
            light.add(phases)
    assert table == {(0, 1): 1, (lower_mass, 1): 2, (high_mass, 1): 2,
                     (upper_mass, -1): 10, (8, -1): 5, (4, 1): 10}
    assert light == {(a,)*5 for a in (0, 6, 24)}
    assert sum(sign*count for (_, sign), count in table.items()) == 0
    return table


def sharp_gap_bounds():
    eps = Q(23, 1000)
    assert 0 < (upper_mass-lower_mass)/2 < eps
    spread = 2*sp.sqrt(2)*eps
    # h^2 bounds the square of every absolute real part from below.
    h2 = Q(3, 10)**2 - spread
    h2_four = Q(3, 8)**2 - spread
    assert h2 > eps and 9*h2_four > 25*eps
    assert Q(25, 2)*eps**2/h2 + eps < 4-upper_mass
    for signs in product((-1, 1), repeat=5):
        d = sum(signs)
        shift = sp.eye(5) - sp.ones(5, 1)*sp.Matrix([signs])/d
        centered = sp.eye(5)-sp.ones(5)/5
        assert shift*centered == shift
        assert max((shift.T*shift).eigenvals()) == Q(25, d*d)
    assert Q(1, 50) < Q(139, 40)-lower_mass < upper_mass-Q(139, 40)


def gw_algebra():
    # gamma5 anticommutes with this Clifford vector of arbitrary length.
    b, s, r = sp.symbols("b s r", real=True)
    gamma5 = sp.diag(1, -1)
    vector = sp.Matrix([[0, s], [s, 0]])
    unit = sp.eye(2)
    x = b*unit + sp.I*vector
    assert gamma5*x*gamma5 == x.H
    assert x.H*x == (b*b+s*s)*unit
    operator = unit + x/r
    residual = (gamma5*operator + operator*gamma5 - operator*gamma5*operator)*r**2
    assert residual.applyfunc(lambda t: sp.expand(t).subs(r*r, b*b+s*s)) == sp.zeros(2)


def finite_locality_bound():
    # For a general degree the coefficients past the midpoint decrease
    # term by term. These examples also check the endpoint convention N=0.
    for degree in (0, 1, 2, 7, 24):
        coeff = [Q(comb(2*k, k), 4**k) for k in range(degree+1)]
        squared = [sum(coeff[i]*coeff[k-i] for i in range(max(0, k-degree),
                                                          min(k, degree)+1))
                   for k in range(2*degree+1)]
        assert squared[:degree+1] == [1]*(degree+1)
        tail = [a-b for a, b in zip(squared[degree:], squared[degree+1:]+[0])]
        assert min(tail) > 0 and sum(tail) == 1


def flavor_symmetry():
    labels = {0, 1, 4}
    affine = [(a, b) for a in range(1, 5) for b in range(5)
              if {(a*x+b) % 5 for x in labels} == labels]
    assert affine == [(1, 0), (4, 0)]


if __name__ == "__main__":
    table = cone_table()
    sharp_gap_bounds()
    gw_algebra()
    finite_locality_bound()
    flavor_symmetry()
    for (mass, sign), count in sorted(table.items(), key=lambda row: row[0][0]):
        print("Wilson value, cones, orientation:", mass, count, sign)
    print("Uniform sharp gap, finite-range approximation and order-five selection verified.")
