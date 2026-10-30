from fractions import Fraction
from itertools import product

I = (0, 1)
ZERO = (Fraction(0), Fraction(0))

def gadd(a, b): return (a[0] + b[0], a[1] + b[1])
def gneg(a): return (-a[0], -a[1])
def gmul(a, b): return (a[0]*b[0] - a[1]*b[1], a[0]*b[1] + a[1]*b[0])
def gscale(q, a): return (q*a[0], q*a[1])

def pprod(a, b):
    # Pauli multiplication in the order I,X,Y,Z, with XY=iZ, YZ=iX, ZX=iY.
    table = {
        ('I','I'): ('I',0), ('I','X'): ('X',0), ('I','Y'): ('Y',0), ('I','Z'): ('Z',0),
        ('X','I'): ('X',0), ('X','X'): ('I',0), ('X','Y'): ('Z',1), ('X','Z'): ('Y',3),
        ('Y','I'): ('Y',0), ('Y','X'): ('Z',3), ('Y','Y'): ('I',0), ('Y','Z'): ('X',1),
        ('Z','I'): ('Z',0), ('Z','X'): ('Y',1), ('Z','Y'): ('X',3), ('Z','Z'): ('I',0),
    }
    out, phase = [], 0
    for x, y in zip(a, b):
        z, k = table[(x, y)]
        out.append(z); phase = (phase + k) % 4
    return tuple(out), ((1,0),(0,1),(-1,0),(0,-1))[phase]

def add(A, B):
    C = dict(A)
    for w, c in B.items():
        C[w] = gadd(C.get(w, ZERO), c)
        if C[w] == ZERO: del C[w]
    return C

def scale(q, A): return {w: gscale(q, c) for w,c in A.items() if gscale(q,c) != ZERO}
def gscale_gaussian(q, A): return {w: gmul(q, c) for w,c in A.items() if gmul(q,c) != ZERO}

def mul(A, B):
    C = {}
    for wa, ca in A.items():
        for wb, cb in B.items():
            w, ph = pprod(wa, wb)
            C[w] = gadd(C.get(w, ZERO), gmul(ca, gmul(ph, cb)))
    return {w:c for w,c in C.items() if c != ZERO}

def comm(A, B): return add(mul(A,B), scale(-1, mul(B,A)))
def anticom(A, B): return add(mul(A,B), mul(B,A))

def one(n): return {('I',)*n: (Fraction(1), Fraction(0))}
def pauli(n, entries):
    w = ['I']*n
    for j, x in entries.items(): w[j-1] = x
    return {tuple(w): (Fraction(1), Fraction(0))}

def ham(m, bs):
    n = m + 2
    H = {}
    for j, b in enumerate(bs, 1):
        H = add(H, scale(Fraction(b), pauli(n, {j:'Z', j+1:'Z', j+2:'X'})))
    return H

def krylov(m, bs):
    H = ham(m, bs)
    eta = [pauli(m+2, {1:'X'})]
    for _ in range(3):
        eta.append(scale(Fraction(1,2), gscale_gaussian(I, comm(H, eta[-1]))))
    return eta

def scalar(A, n):
    assert len(A) <= 1
    if not A: return Fraction(0)
    w, c = next(iter(A.items()))
    assert w == ('I',)*n and c[1] == 0
    return c[0]

def gram(m, bs):
    E = krylov(m, bs)
    n = m+2
    return [[scalar(anticom(E[j], E[k]), n) for k in range(4)] for j in range(4)]

def delta(C):
    return C[1][1]**2 * (C[0][0]*C[2][2] - C[1][1]**2) - C[0][0]**2 * (C[1][1]*C[3][3] - C[2][2]**2)

def det2(a,b,c,d): return a*d-b*c

def check(bs, expected, label):
    C = gram(len(bs), bs)
    got = delta(C)
    assert got == expected, (label, got, expected, C)
    return C

# Fidelity: the paper's six unit couplings give the displayed block and Δ=-16.
C6 = check([1,1,1,1,1,1], Fraction(-16), 'six-unit control')
assert C6 == [[2,0,-2,0],[0,2,0,-6],[-2,0,6,0],[0,-6,0,24]]

# The proposed family has Δ=0 with every coupling positive.
for m in [5,6,7,10,15]:
    bs = [3,1,1,1,4] + [2+i for i in range(m-5)]
    C = check(bs, Fraction(0), f'family m={m}')
    x, R, S = bs[0]**2, bs[1]**2 + bs[2]**2, bs[0]**2 + bs[1]**2 + bs[2]**2
    T = bs[3]**2*R + bs[2]**2*bs[4]**2
    assert (bs[0]**2-bs[3]**2)*R - bs[2]**2*bs[4]**2 == 0
    assert det2(C[0][0], C[0][2], C[2][0], C[2][2]) > 0
    assert det2(C[1][1], C[1][3], C[3][1], C[3][3]) > 0

print('six-unit C =', C6)
print('family witness m=7, b=(3,1,1,1,4,1,1), Δ=0')
print('ALL_CHECKS_OK')
