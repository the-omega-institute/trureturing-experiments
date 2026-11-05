#!/usr/bin/env python3
"""Exact checks for the operational conjecture following SLT Theorem 2.

No numerical approximation is used. Symbolic identities check the parameter family;
an independent matrix-eigenvalue calculation checks its rational witness. This is
not a Lean proof, an SDP feasibility claim, or a key-rate calculation.
"""
import sympy as S
from itertools import product

Q = S.Rational
I = S.eye(2)
X = S.Matrix([[0, 1], [1, 0]])
Z = S.diag(1, -1)
kron = S.kronecker_product


def eq(x, y):
    if isinstance(x, S.MatrixBase):
        assert x.shape == y.shape
        assert all(S.simplify(z) == 0 for z in x-y), (x, y)
    else:
        assert S.simplify(x-y) == 0, (x, y)


def ket(i):
    return I[:, i]


def proj(v):
    return v*v.conjugate().T


def flag(P, Qm):
    return (kron(P, proj(ket(0)))+kron(Qm, proj(ket(1))))/2


def distance(A, B):
    H = A-B
    eq(H, H.conjugate().T)
    return S.simplify(sum(S.Abs(x)*m for x, m in H.eigenvals().items())/2)


# Generic family: a^2+b^2=1, s^2+t^2=1, c^2+r^2=1.
a, b, s, t, c, r = S.symbols('a b s t c r', real=True)
e0 = S.Matrix([1, 0])
e1 = S.Matrix([s, t])
R = S.Matrix([[c, -r], [r, c]])
psi = S.Matrix.vstack(a*e0, b*R*e0, -b*R*e1, a*e1)/S.sqrt(2)
relations = S.groebner([a*a+b*b-1, t*t+s*s-1, c*c+r*r-1], b, t, r, a, s, c, domain=S.QQ)


def reduced(x):
    return S.factor(relations.reduce(S.expand(x))[1])


assert reduced((psi.T*psi)[0]-1) == 0
U = S.Matrix([[a, -b], [b, a]])
# Multiply the CHSH observable by sqrt(1+s^2) before polynomial reduction.
Bplus_scaled = U*(Z+s*X)*U.T
Bminus_scaled = U*(Z-s*X)*U.T
chsh_scaled = (psi.T*kron(kron(Z, Bplus_scaled+Bminus_scaled)
                                   + kron(X, Bplus_scaled-Bminus_scaled), I)*psi)[0]
expected_scaled = 2*(1+s*s)*((a*a-b*b)**2+4*a*a*b*b*c)
assert reduced(chsh_scaled-expected_scaled) == 0
for v in [e0, e1]:
    P = proj(v)
    assert reduced((P-proj(R*v)).trace()) == 0
    assert reduced(S.trace((P-proj(R*v))**2)-2*r*r) == 0
print('FAMILY: normalization, pure-state squared distances, and CHSH identity checked symbolically')

# Strict nondegenerate rational point. Tensor order is A,B,E; flag order E,T.
a, b, s, t = Q(4, 5), Q(3, 5), Q(3, 5), Q(4, 5)
c, r = Q(9999, 10001), Q(200, 10001)
e0, e1 = S.Matrix([1, 0]), S.Matrix([s, t])
R = S.Matrix([[c, -r], [r, c]])
eq(R.T*R, I)
psi = S.Matrix.vstack(a*e0, b*R*e0, -b*R*e1, a*e1)/S.sqrt(2)
eq((psi.T*psi)[0], 1)
P0, P1 = proj(e0), proj(e1)
tau0, tau1 = flag(P0, P1), flag(P1, P0)
V = kron(R, I)
eta0, eta1 = V*tau0*V.T, V*tau1*V.T
states = {(0, 0): tau0, (1, 1): tau1, (0, 1): eta0, (1, 0): eta1}
weights = {}
raw = {}
for u, v in product(range(2), repeat=2):
    vec = psi[2*(2*u+v):2*(2*u+v)+2, 0]
    p = S.simplify((vec.T*vec)[0])
    assert p > 0
    weights[u, v], raw[u, v] = p, proj(vec)/p
    eq(p, (a*a if u == v else b*b)/2)
for u, v in product(range(2), repeat=2):
    eq(weights[u, v], weights[1-u, 1-v])
    eq(states[u, v], flag(raw[u, v], raw[1-u, 1-v]))
    eq(S.trace(states[u, v]), 1)
    assert all(x >= 0 for x in states[u, v].eigenvals())
assert all(states[x] != states[y] for x, y in product(states, repeat=2) if x != y)
print('WITNESS: all outcome probabilities positive; Eq. (6) flags match; all four density states distinct')

# Fidelity to the source: its Lemma 1 is a flag-swap unitary equivalence.
# For every length-1, length-2, and length-3 message, check the actual matrices.
flag_flip = kron(I, X)
eq(flag_flip*tau0*flag_flip.T, tau1)
for k in range(1, 4):
    all0, all1 = kron(*([tau0]*k)), kron(*([tau1]*k))
    for msg in product(range(2), repeat=k):
        W = kron(*[flag_flip if z else S.eye(4) for z in msg])
        rho0 = kron(*[tau1 if z else tau0 for z in msg])
        rho1 = kron(*[tau0 if z else tau1 for z in msg])
        eq(W*rho0*W.T, all0)
        eq(W*rho1*W.T, all1)
print('SOURCE FIDELITY: Lemma 1 unitary equivalence reproduced for every message of length 1, 2, 3')

# Exact spectral calculations, independent of the symbolic squared-distance check.
dcorrect = distance(tau0, tau1)
d0, d1 = distance(tau0, eta0), distance(tau1, eta1)
eq(dcorrect, Q(4, 5))
eq(d0, r)
eq(d1, r)
obstruction = S.simplify(dcorrect-d0-d1)
eq(obstruction, Q(38004, 50005))
assert obstruction > 0
print(f'TRACE DISTANCES: D(correct0,correct1)={dcorrect}; D(correct_i,error_i)={d0}')
print(f'CPTP OBSTRUCTION: total correct-state disturbance + residual error distance >= {obstruction} > 0')

# Rank-one flagged states have spectrum {1/2,1/2,0,0}; all positive powers
# are 2^(-u) times their support projector. Thus Q(u) is constant on 0<u<1.
H0, H1 = 2*tau0, 2*tau1
eq(H0*H0, H0)
eq(H1*H1, H1)
u = S.symbols('u', real=True)
q_at_u = S.simplify(2**(-u)*2**(-(1-u))*S.trace(H0*H1))
eq(q_at_u, s*s)
# Root fidelity is Tr sqrt(sqrt(tau0)*tau1*sqrt(tau0)).
sqrt_tau0 = S.sqrt(2)*tau0
eq(sqrt_tau0*sqrt_tau0, tau0)
sandwich = S.simplify(sqrt_tau0*tau1*sqrt_tau0)
fidelity = S.simplify(sum(S.sqrt(x)*m for x,m in sandwich.eigenvals().items()))
eq(fidelity, Q(3, 5))
assert fidelity**2 <= q_at_u <= fidelity  # source's two stated comparisons
error = weights[0, 1]+weights[1, 0]
threshold = S.simplify(error/(1-error))
eq(error, Q(9, 25))
eq(threshold, Q(9, 16))
assert q_at_u < threshold
print(f'SOURCE FIDELITY: F={fidelity}, Q={q_at_u}; F^2 <= Q <= F reproduced')
print(f'CHERNOFF REGION: Q={q_at_u} < epsilon/(1-epsilon)={threshold}, epsilon={error}')

U = S.Matrix([[a, -b], [b, a]])
B1, B2 = U*(Z+s*X)*U.T/S.sqrt(1+s*s), U*(Z-s*X)*U.T/S.sqrt(1+s*s)
for O in [Z, X, B1, B2]:
    eq(O.T, O)
    eq(O*O, I)
chsh = S.simplify((psi.T*kron(kron(Z, B1+B2)+kron(X, B1-B2), I)*psi)[0])
expected = 2*S.sqrt(34)/5*(1-Q(1152, 6250625))
eq(chsh, expected)
assert chsh > 0 and S.simplify(chsh**2-4) > 0
print(f'BELL VIOLATION: CHSH={chsh}; CHSH^2-4={S.factor(chsh**2-4)} > 0')
print('ALL_CHECKS_OK')
