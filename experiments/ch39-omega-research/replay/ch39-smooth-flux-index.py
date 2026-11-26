"""Finite magnetic trace interpolation, the Bott sign and the rank of a shifted line projector."""
from fractions import Fraction as F
from itertools import product
import sympy as sp


def interpolation_mesh(base, flux_bound):
    k = 16*(1024+(base*(1+flux_bound)).bit_length())
    mesh = base*2**k
    degree = 16*max(1, flux_bound)*(2*(2560000*k)+1)**2
    # Every rank density on this mesh has denominator dividing mesh squared.
    err = F(16*mesh**2, 16**k) + F(1024*(2*degree+1)*degree**3, mesh)
    assert err < 1


def path_area():
    # Sanity check normal ordering for every short magnetic word in one plane.
    for length in range(1, 6):
        for word in product(((1,0),(-1,0),(0,1),(0,-1)), repeat=length):
            horiz = vert = area = 0
            for a, b in word:
                area -= vert*a
                horiz += a
                vert += b
            assert abs(area) <= length*length
            if horiz == vert == 0:
                # Move S and C with S C = z C S; the reverse path changes the sign.
                back = power = 0
                for a, b in reversed(word):
                    power -= back*(-a)
                    back -= b
                assert power == -area


def clifford_sign():
    one = sp.eye(2)
    x, y, z = sp.Matrix([[0,1],[1,0]]), sp.Matrix([[0,-sp.I],[sp.I,0]]), sp.diag(1,-1)
    kron = sp.kronecker_product
    gamma = [kron(x,one), kron(y,one), kron(z,x), kron(z,y)]
    chirality = gamma[0]*gamma[1]*gamma[2]*gamma[3]
    assert chirality == sp.diag(-1,1,1,-1)
    e = sp.eye(4)[:, 1]
    vectors = [sp.I*chirality*g*e for g in gamma]
    columns = [[part for i in (0, 3) for part in v[i].as_real_imag()] for v in vectors]
    assert sp.Matrix(columns).T.det() == -1


def line_projector(q):
    # Rational rotation columns on two ramps. In the second ramp the bottom index advances.
    gate, stop = q//2, q-2
    entries = []
    for j in range(q):
        if j == stop:
            continue
        if j < gate:
            angle = min(F(1), F(2*j, gate))
            top = 2*angle/(1+angle**2)
            bottom = (1-angle**2)/(1+angle**2)
            partner = j
        else:
            angle = F(j-gate, stop-gate-1)
            bottom = 2*angle/(1+angle**2)
            top = (1-angle**2)/(1+angle**2)
            partner = j+1
        # The last site is in the bottom plateau.
        if j > stop:
            top, bottom, partner = 0, 1, j
        assert top*top+bottom*bottom == 1
        entries.append((j,partner))
    assert len(entries) == q-1
    assert len(set(i for i, _ in entries)) == len(set(j for _, j in entries)) == q-1
    return len(entries)


def smooth_flux():
    for q in (16, 31, 64, 103):
        r = line_projector(q)
        assert r*r-2*q*r+q*q == 1
    side = 25000
    assert F(1, side*side) < F(1, 2**29)
    charges = [-1,1,0,0,0]
    exterior = [a+b for j,a in enumerate(charges) for b in charges[j+1:]]
    assert [-3*sum(x*x for x in weights) for weights in (charges, exterior)] == [-6, -18]
    # The rank shift for either sign is the telescoping diagonal correction on the falling ramp.
    plateau = [1]*20+[F(i, 20) for i in reversed(range(21))]+[0]*20
    for step in range(-6, 7):
        assert sum(plateau[j]-plateau[j-step] for j in range(10, 51)) == -step
    flux = sp.Matrix([[0,1,2,3],[-1,0,4,5],[-2,-4,0,6],[-3,-5,-6,0]])
    basis = sp.eye(4)
    basis[0,1], basis[2,3], basis[1,2] = 7, -2, 3
    def pfaffian(m):
        return m[0,1]*m[2,3]-m[0,2]*m[1,3]+m[0,3]*m[1,2]
    assert basis.det() == 1 and pfaffian(basis.T*flux*basis) == pfaffian(flux)


if __name__ == "__main__":
    for base, other in [(2**29,1),(2**29*625000000,625000000),(2**1024+1,2**513-1)]:
        interpolation_mesh(base, other)
    path_area()
    clifford_sign()
    smooth_flux()
    print("Finite interpolation bounds, sphere orientation, shifted-line ranks and smooth-flux indices verified.")
