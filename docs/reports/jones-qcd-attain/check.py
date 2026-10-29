#!/usr/bin/env python3
"""Exact checks for Jones arXiv:2610.16270v2, Section 13(1).

No floating point, numerical eigensolver, or tolerance is used. These finite
checks do not prove the universal nonattainment statement.
"""
import hashlib
from pathlib import Path
import sympy as s

Q = s.Rational
I = s.eye(2)
Z = s.zeros(2)

def eq(a, b):
    if isinstance(a, s.MatrixBase):
        assert a.shape == b.shape
        assert all(s.simplify(x) == 0 for x in a - b), (a, b)
    else:
        assert s.simplify(a - b) == 0, (a, b)

def nonneg(x):
    x = s.simplify(x)
    assert x.is_nonnegative is True, x

def psd2(A):
    eq(A, A.H)
    nonneg(A[0, 0]); nonneg(A[1, 1]); nonneg(A.det())

def norm2(A):
    """Exact Hilbert operator norm of a real symmetric 2 by 2 matrix."""
    eq(A, A.T)
    gap = s.sqrt((A[0, 0] - A[1, 1])**2 + 4*A[0, 1]**2)
    return (s.Abs(s.trace(A)) + gap)/2

def choi(phi):
    J = s.zeros(4)
    for i in range(2):
        for j in range(2):
            E = s.zeros(2); E[i, j] = 1
            B = phi(E)
            for y in range(2):
                for z in range(2):
                    J[2*y+i, 2*z+j] = B[y, z]
    return J

def mean_diagonal(A, B):
    assert A.is_diagonal() and B.is_diagonal()
    for x in list(A.diagonal()) + list(B.diagonal()): nonneg(x)
    return s.diag(*(s.sqrt(A[i, i]*B[i, i]) for i in range(A.rows)))

def partial_trace_output(J):
    return s.Matrix(2, 2, lambda i, j: sum(J[2*y+i, 2*y+j] for y in range(2)))

def replacement(D):
    return lambda X: s.trace(X)*D

def readout(D, E):
    H = partial_trace_output(mean_diagonal(choi(replacement(D)), choi(replacement(E))))
    assert H.is_diagonal()
    eq(H[0, 0], H[1, 1])
    return s.simplify(H[0, 0])

print('Arithmetic: SymPy exact integers, rationals, and algebraic square roots only')
# Source fidelity: Proposition 3.2 (label impossibility proposition), choosing
# orthogonal pure states. The split pair has fidelity zero; the equal split one.
P = s.diag(1, 0); R = s.diag(0, 1)
eq(readout(P, R), 0)
eq(readout(I/2, I/2), 1)
print('SOURCE_FIDELITY: Proposition 3.2: Fhat(Tr(.)P,Tr(.)R)=0; Fhat(Omega/2,Omega/2)=1')

# Four nonzero CPTP replacement channels. Compute each unnormalised Choi
# matrix directly from its action on matrix units, not from a supplied table.
parameters = [Q(1, 16), Q(1, 8), Q(1, 4), Q(1, 3), Q(1, 2)]
for d in parameters:
    D = s.diag(d, 1-d); E = I-D
    for output in [D, E, R, P]:
        eq(s.trace(output), 1)
        assert output != Z
        kraus = []
        for y in range(2):
            for x in range(2):
                K = s.zeros(2); K[y, x] = s.sqrt(output[y, y]); kraus.append(K)
        eq(sum((K.H*K for K in kraus), Z), I)
        for i in range(2):
            for j in range(2):
                X = s.zeros(2); X[i, j] = 1
                eq(sum((K*X*K.H for K in kraus), Z), replacement(output)(X))
    eq(choi(replacement(D))+choi(replacement(E)), s.eye(4))
    eq(choi(replacement(R))+choi(replacement(P)), s.eye(4))
    # Matrix-unit Kraus basis has W=I_4 in (output,input) order.
    T1 = choi(replacement(D)); T2 = choi(replacement(R))
    eq(T1-T2, s.diag(d, d, -d, -d))
    kappa = 2; t = d
    f1, f2 = readout(D, E), readout(R, P)
    eq(f1, 2*s.sqrt(d*(1-d))); eq(f2, 0)
    eq(f1, kappa*s.sqrt(t*(1-t)))
    nonneg(kappa*s.sqrt(t)-f1)
    assert s.simplify(kappa*s.sqrt(t)-f1).is_positive is True
    eq((f1/(kappa*s.sqrt(t)))**2, 1-d)
    print(f'CHANNEL: delta={d}; kappa=2; t={d}; Fhat1={f1}; Fhat2=0; old_bound={2*s.sqrt(d)}; STRICT')

# Generic symbolic anticommutator identity, without any commutativity assumption.
E, F = s.symbols('E F', commutative=False)
A, B = E-F, 1-E-F
eq(s.expand((A*B+B*A)/2-((E-E**2)-(F-F**2))), 0)
print('IDENTITY: q(E)-q(F)=(AB+BA)/2, checked in the free noncommutative algebra')

# All pairs from a fixed collection of noncommuting rational effects. Spectral
# data are specified exactly; their PSD and effect conditions are also checked.
rotations = [I, s.Matrix([[Q(3,5),-Q(4,5)],[Q(4,5),Q(3,5)]])]
spectra = [(0,0),(0,Q(1,4)),(0,Q(1,2)),(0,1),(Q(1,4),Q(3,4)),(Q(1,2),1),(1,1)]
effects=[]
for U in rotations:
    eq(U.T*U,I)
    for a,b in spectra:
        E=U*s.diag(a,b)*U.T
        g=U*s.diag(s.sqrt(a*(1-a)),s.sqrt(b*(1-b)))*U.T
        psd2(E); psd2(I-E); eq(g*g,E-E*E)
        effects.append((E,g))
count=noncommuting=0
for E,gE in effects:
    for F,gF in effects:
        t=s.simplify(norm2(E-F)); nonneg(t); nonneg(1-t)
        small=s.simplify(Q(1,2)-t).is_nonnegative
        assert small is not None
        d=s.simplify(t*(1-t)) if small else Q(1,4)
        D=E-E*E-F+F*F
        psd2(d*I-D); psd2(d*I+D)
        K=gE-gF
        psd2(d*I-K*K)
        if t != 0:
            assert s.simplify(t-d).is_positive is True
        if E*F != F*E: noncommuting+=1
        count+=1
print(f'EFFECTS: {count} exact ordered pairs; {noncommuting} noncommuting; quadratic and square-root bounds pass')
print('LIMIT: channel ratio squared = 1-delta, whose limit at delta down to 0 is 1')
print('Scope: finite exact checks and symbolic identities; no universal Lean proof')
print('check.py sha256='+hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
print('ALL_CHECKS_OK')
