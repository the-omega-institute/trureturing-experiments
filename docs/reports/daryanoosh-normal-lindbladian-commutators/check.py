#!/usr/bin/env python3
"""Exact checks for arXiv:2608.04775v1, Sections VIII and IX.

Source fidelity: Section VIII's driven infinite-temperature thermal qubit.
Candidate: a symbolic five-parameter GKSL family with three nonzero sectors.
This is algebraic evidence, not Lean kernel verification or a novelty claim.
Requires Python 3 and SymPy; uses no floating-point evidence.
"""
import sys
import sympy as s

I = s.I
ONE = s.eye(2)
X = s.Matrix([[0, 1], [1, 0]])
Y = s.Matrix([[0, -I], [I, 0]])
Z = s.diag(1, -1)
PAULI = (ONE, X, Y, Z)


def zero(A):
    return all(s.expand(v) == 0 for v in A)


def equal(A, B, label):
    if not zero(A - B):
        raise AssertionError(f'{label}: {s.simplify(A - B)}')


def comm(A, B):
    return A * B - B * A


def parts(H, jumps):
    n = H.rows
    loss = sum((L.H * L for L in jumps), s.zeros(n))
    effective = H - I * loss / 2

    def smooth(rho):
        return -I * (effective * rho - rho * effective.H)

    def jump(rho):
        return sum((L * rho * L.H for L in jumps), s.zeros(n))

    def gksl(rho):
        return -I * comm(H, rho) + sum(
            (L * rho * L.H - (L.H * L * rho + rho * L.H * L) / 2
             for L in jumps), s.zeros(n))
    return loss, smooth, jump, gksl


def matrix_units(n):
    ans = []
    for i in range(n):
        for j in range(n):
            E = s.zeros(n)
            E[i, j] = 1
            ans.append(E)
    return ans


def units_matrix(F, n=2):
    # E_ij in row-major order are orthonormal for Tr(A† B).
    return s.Matrix.hstack(*(s.Matrix(list(F(E))) for E in matrix_units(n)))


def pauli_matrix(F):
    # All four basis elements have norm sqrt(2), so this is also the
    # matrix in the orthonormal basis (I,X,Y,Z)/sqrt(2).
    return s.Matrix(4, 4, lambda i, j: s.expand(s.trace(PAULI[i].H * F(PAULI[j])) / 2))


def sectors(S, J):
    return tuple(A.applyfunc(s.expand) for A in (
        comm(S, S.H), comm(J, J.H), comm(S, J.H) + comm(J, S.H)))


def fidelity():
    q = s.symbols('q', positive=True)
    nu = s.symbols('nu', real=True)
    plus, minus = (X + I * Y) / 2, (X - I * Y) / 2
    H = nu * Z / 2
    R, smooth, jump, full = parts(H, [q * plus, q * minus])
    expected = s.Matrix([[0, 0, 0, 0], [0, -q**2, -nu, 0],
                         [0, nu, -q**2, 0], [0, 0, 0, -2*q**2]])
    equal(R, q**2 * ONE, 'source thermal uniform loss')
    equal(pauli_matrix(full), expected, 'source thermal Pauli matrix')
    equal(comm(expected, expected.H), s.zeros(4), 'source thermal normality')
    for D in sectors(pauli_matrix(smooth), pauli_matrix(jump)):
        equal(D, s.zeros(4), 'source thermal sector is zero')
    print('FIDELITY: driven thermal qubit matrix, normality, and all three zero sectors OK (gamma_th=q^2)')
    print('thermal_matrix =', expected.tolist())


def candidate():
    r, t, a, b, h = s.symbols('r t a b h', positive=True)
    jumps = [r * X + a * ONE, r * Y + I * b * ONE, t * Z]
    H = h * Z - r * b * Y
    R, smooth, jump, full = parts(H, jumps)
    equal(H.H, H, 'candidate Hamiltonian Hermitian')
    c = a**2 + b**2 + 2*r**2 + t**2
    equal(R, c * ONE + 2*a*r*X, 'candidate loss')
    S, J, G = map(pauli_matrix, (smooth, jump, full))
    expected = s.Matrix([[0, 0, 0, 0],
                         [0, -2*(r**2+t**2), -2*h, 0],
                         [0, 2*h, -2*(r**2+t**2), 0],
                         [0, 0, 0, -4*r**2]])
    equal(G, expected, 'candidate GKSL reduction')
    equal(G, S + J, 'candidate decomposition')
    equal(comm(G, G.H), s.zeros(4), 'candidate normality')
    ds, dj, mixed = sectors(S, J)
    U, V, W = s.zeros(4), s.zeros(4), s.zeros(4)
    U[0, 2] = U[2, 0] = 1
    V[0, 3] = V[3, 0] = 1
    W[1, 3] = W[3, 1] = 1
    equal(ds, -8*a*r*h*U - 8*a*b*r**2*V, 'candidate D_S')
    equal(dj, -8*a*b*r**2*V + 8*b*r*(t**2-r**2)*W, 'candidate D_J')
    equal(mixed, 8*a*r*h*U + 16*a*b*r**2*V - 8*b*r*(t**2-r**2)*W, 'candidate D_SJ')
    equal(ds + dj + mixed, s.zeros(4), 'three-sector cancellation')
    witnesses = (s.factor(ds[0, 2]), s.factor(dj[0, 3]), s.factor(mixed[0, 3]))
    if not (witnesses[0].is_negative and witnesses[1].is_negative and witnesses[2].is_positive):
        raise AssertionError('nonvanishing witness sign not certified')

    # Independent coordinate check: act on every E_ij, then independently
    # assemble the superoperators via Kronecker left/right multiplication.
    Sc, Jc, Gc = map(units_matrix, (smooth, jump, full))
    equal(Gc, Sc + Jc, 'matrix-unit decomposition')
    equal(comm(Gc, Gc.H), s.zeros(4), 'matrix-unit normality')
    Heff = H - I*R/2
    Sk = -I * (s.kronecker_product(Heff, ONE) - s.kronecker_product(ONE, Heff.conjugate()))
    Jk = sum((s.kronecker_product(L, L.conjugate()) for L in jumps), s.zeros(4))
    equal(Sc, Sk, 'independent Kronecker smooth matrix')
    equal(Jc, Jk, 'independent Kronecker jump matrix')
    P = s.Matrix.hstack(*(s.Matrix(list(B)) for B in PAULI))
    equal(P.H * P, 2*s.eye(4), 'Pauli Gram matrix')
    for Dc, Dp in zip(sectors(Sc, Jc), (ds, dj, mixed)):
        equal(P.H * Dc * P / 2, Dp, 'sector coordinate bridge')
        if zero(Dc):
            raise AssertionError('a matrix-unit sector vanishes')

    lam = s.symbols('lambda')
    characteristic = lam*(lam+4*r**2)*((lam+2*(r**2+t**2))**2+4*h**2)
    if s.expand(G.charpoly(lam).as_expr()-characteristic) != 0:
        raise AssertionError('characteristic polynomial mismatch')
    for E in matrix_units(2):
        equal(full(E), -I*h*comm(Z,E)+r**2*(X*E*X-E)+r**2*(Y*E*Y-E)+t**2*(Z*E*Z-E), 'operator reduction on E_ij')
    equal(full(ONE/2), s.zeros(2), 'stationary density I/2')
    # Exclude hidden inexact arithmetic anywhere in the evidence matrices.
    if any(A.has(s.Float) for A in [H,R,S,J,G,Sc,Jc,Gc,ds,dj,mixed]):
        raise AssertionError('floating point expression present')
    print('FAMILY: exact polynomial identities valid for positive r,t,a,b,h; requested subdomain is 0<r<t')
    print('full_generator =', G.tolist())
    print('D_S =', ds.tolist())
    print('D_J =', dj.tolist())
    print('D_SJ =', mixed.tolist())
    print('nonzero_entries =', witnesses)
    print('COMPUTATIONAL_BASIS: normality, decomposition, nonzero sectors and independent Kronecker assembly OK')
    print('CHARACTERISTIC_POLYNOMIAL:', s.factor(characteristic))
    print('STATIONARY: I/2 OK')


if __name__ == '__main__':
    print('Python', sys.version.split()[0], 'SymPy', s.__version__)
    fidelity()
    candidate()
    print('ALL_CHECKS_OK')
