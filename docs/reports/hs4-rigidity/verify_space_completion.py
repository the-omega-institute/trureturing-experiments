#!/usr/bin/env python3
"""Finite exact checks and seeded tests for HS4 spatial completion.
No Lean replay, universal numerical certification, or physical realization.
"""
from __future__ import annotations
import itertools
import json
import platform
import numpy as np
import sympy as s

exact = numerical = negative = 0

def eq(a, b, label):
    global exact
    difference = a - b
    if isinstance(difference, s.MatrixBase):
        ok = all(s.simplify(x) == 0 for x in difference)
    else:
        ok = s.simplify(difference) == 0
    if not ok:
        raise RuntimeError(label)
    exact += 1

def near(a, b, label, tol=2e-11):
    global numerical
    if np.max(np.abs(np.asarray(a) - np.asarray(b))) > tol:
        raise RuntimeError(label)
    numerical += 1

def le(a, b, label, tol=2e-11):
    global numerical
    if a > b + tol:
        raise RuntimeError(label)
    numerical += 1

I = s.eye(2)
X = s.Matrix([[0, 1], [1, 0]])
Y = s.Matrix([[0, -s.I], [s.I, 0]])
Z = s.diag(1, -1)
paulis = [X, Y, Z]

def kron(*args):
    return s.kronecker_product(*args)

def local(op, k):
    return kron(*[op if i == k else I for i in range(4)])

def swap(i, j):
    m = s.zeros(16)
    for n in range(16):
        bits = [(n >> (3-k)) & 1 for k in range(4)]
        bits[i], bits[j] = bits[j], bits[i]
        out = sum(b << (3-k) for k, b in enumerate(bits))
        m[out, n] = 1
    return m

def reduce_pair(rho, keep):
    other = [k for k in range(4) if k not in keep]
    m = s.zeros(4)
    for a, b, t in itertools.product(range(4), range(4), range(4)):
        p, q = [0]*4, [0]*4
        for k, site in enumerate(keep):
            p[site], q[site] = (a >> (1-k)) & 1, (b >> (1-k)) & 1
        for k, site in enumerate(other):
            p[site] = q[site] = (t >> (1-k)) & 1
        pi = sum(v << (3-k) for k, v in enumerate(p))
        qi = sum(v << (3-k) for k, v in enumerate(q))
        m[a, b] += rho[pi, qi]
    return m

e0, e1 = s.zeros(16, 1), s.zeros(16, 1)
for i in [3, 12]:
    e0[i] = 1/s.sqrt(3)
for i in [5, 6, 9, 10]:
    e0[i] = -1/(2*s.sqrt(3))
for i, sign in [(5, 1), (6, -1), (9, -1), (10, 1)]:
    e1[i] = s.Rational(sign, 2)
E = e0.row_join(e1)
eq(E.H*E, I, 'isometry')
for p in paulis:
    eq(sum((local(p, k) for k in range(4)), s.zeros(16))*E, s.zeros(16, 2), 'closure')
logical_swaps = [Z, -Z/2-s.sqrt(3)*X/2, -Z/2+s.sqrt(3)*X/2]
pairs = list(itertools.combinations(range(4), 2))
which = {(0,1):0, (0,2):1, (0,3):2, (1,2):2, (1,3):1, (2,3):0}
u, v, w = s.symbols('u v w', real=True)
logical = (I+u*X+v*Y+w*Z)/2
physical = E*logical*E.H
singlet = s.Matrix([0, 1, -1, 0])/s.sqrt(2)
ps = singlet*singlet.H
values = [w, -w/2-s.sqrt(3)*u/2, -w/2+s.sqrt(3)*u/2]
for pair in pairs:
    j = which[pair]
    eq(E.H*swap(*pair)*E, logical_swaps[j], 'swap compression')
    a = values[j]
    werner = (1-a)*ps/2 + (1+a)*(s.eye(4)-ps)/6
    eq(reduce_pair(physical, pair), werner, 'entire pair matrix')
q = s.zeros(16)
for a, b, c in itertools.permutations(range(3)):
    q += s.LeviCivita(a,b,c)*local(paulis[a],0)*local(paulis[b],1)*local(paulis[c],2)/8
qhat = E.H*q*E
eq(qhat, -s.sqrt(3)*Y/4, 'oriented volume')
eq(q*q*E, s.Rational(3,16)*E, 'absolute volume degeneracy')
eq(logical.det(), (1-u*u-v*v-w*w)/4, 'positivity ball')
eq(sum(values), 0, 'boundary plane')
eq(sum(a*a for a in values), 3*(u*u+w*w)/2, 'boundary disk')
eq(logical_swaps[0]*logical_swaps[1]-logical_swaps[1]*logical_swaps[0], -s.I*s.sqrt(3)*Y, 'shape commutator')
f = s.symbols('f', real=True)
gap = ((1+f)*s.log(1+f)+(1-f)*s.log(1-f))/2
eq(s.diff(gap, f, 2), 1/(1-f*f), 'entropy curvature')

eta, kappa, t = s.symbols('eta kappa t', real=True)
def channel(m):
    return (s.trace(m)*I+eta*s.trace(X*m)*X+kappa*s.trace(Y*m)*Y+eta*s.trace(Z*m)*Z)/2
choi = s.zeros(4)
for a,b in itertools.product(range(2),repeat=2):
    unit = s.zeros(2); unit[a,b]=1
    choi += kron(unit,channel(unit))/2
weights = [(1+2*eta+kappa)/4,(1-kappa)/4,(1-2*eta+kappa)/4,(1-kappa)/4]
characteristic = choi.charpoly()
eq(characteristic.as_expr(), s.prod(characteristic.gen-p for p in weights), 'complete Choi spectrum')
D = choi.subs({eta:1,kappa:0})
phi_minus = s.Matrix([0,1,-1,0])/s.sqrt(2)
eq((phi_minus.H*D*phi_minus)[0], -s.Rational(1,4), 'negative Choi witness')
for m in [I,X,Y,Z]:
    eq(channel(m).subs({eta:s.Rational(1,2),kappa:0}), m/2+X*m*X/4+Z*m*Z/4, 'attaining noisy channel')
effects = [(I+X)/4,(I-X)/4,(I+Z)/4,(I-Z)/4]
eq(sum(effects,s.zeros(2)), I, 'classical POVM normalization')
for effect, p in zip(effects,[(1+u)/4,(1-u)/4,(1+w)/4,(1-w)/4]):
    eq(s.trace(logical*effect),p,'classical shape readout')
plus = s.Matrix([1,s.I])/s.sqrt(2); minus=s.conjugate(plus)
ghz = (kron(plus,plus)+kron(minus,minus))/s.sqrt(2)
ghzrho=ghz*ghz.H
for site in [0,1]:
    red=s.zeros(2)
    for a,b,k in itertools.product(range(2),repeat=3):
        red[a,b]+=ghzrho[2*a+k,2*b+k] if site==0 else ghzrho[2*k+a,2*k+b]
    eq(red,I/2,'two-node coherent cat marginal')

false_claims = [qhat == s.zeros(2),
                s.Rational(-1,4) >= 0,
                (I/2)*(I/2) == I/2]
for claim in false_claims:
    if bool(claim):
        raise RuntimeError('negative control unexpectedly accepted')
    negative += 1

rng=np.random.default_rng(20261014)
EN=np.array(E.evalf(),dtype=complex)
PN=np.array(ps.evalf(),dtype=complex)
LN=[np.array(m.evalf(),dtype=complex) for m in logical_swaps]
XN,YN,ZN=[np.array(m,dtype=complex) for m in paulis]
FSTAR=1+np.log2(3)/2

def entropy(p):
    p=np.asarray(p,dtype=float); p=p[p>1e-15]
    return -np.sum(p*np.log2(p))

def local_gap(b):
    vals=np.array([b[2],-b[2]/2-np.sqrt(3)*b[0]/2,-b[2]/2+np.sqrt(3)*b[0]/2])
    fval=np.mean([entropy([(1-a)/2]+[(1+a)/6]*3) for a in vals])
    delta=FSTAR-fval
    r2=b[0]**2+b[2]**2
    le(r2/(4*np.log(2)),delta,'lower deficit bound')
    le(delta,r2/2,'upper deficit bound')
    return delta

for pure in [False,True]:
    for _ in range(200):
        b=rng.normal(size=3); b/=np.linalg.norm(b)
        if not pure: b*=rng.uniform()
        delta=local_gap(b)
        eps=(1-abs(b[1]))/2
        impurity=(1-np.dot(b,b))/2
        le(eps,impurity+2*np.log(2)*delta,'local purity plus entropy budget')
        if pure:
            dist=np.sqrt(2-2*np.sqrt((1+abs(b[1]))/2))
            le(dist,2*np.sqrt(max(0,np.log(2)*delta)),'square root stability')
        rho=EN@((np.eye(2)+b[0]*XN+b[1]*YN+b[2]*ZN)/2)@EN.conj().T
        for pair in pairs:
            order=list(pair)+[k for k in range(4) if k not in pair]
            a=rho.reshape([2]*8).transpose(order+[k+4 for k in order]).reshape(4,4,4,4)
            red=np.einsum('atbt->ab',a)
            j=which[pair]
            sval=np.trace(((np.eye(2)+b[0]*XN+b[1]*YN+b[2]*ZN)/2)@LN[j]).real
            target=(1-sval)*PN/2+(1+sval)*(np.eye(4)-PN)/6
            near(red,target,'sampled pair formula')

nontrivial_budgets = 0
for j in range(200):
    if j < 100:
        a=rng.normal(size=(8,4))+1j*rng.normal(size=(8,4))
        rho=a@a.conj().T;rho/=np.trace(rho)
    else:
        node=np.array(plus,dtype=complex).reshape(2)
        ref=np.kron(np.kron(node,node),node)
        perturb=rng.normal(size=8)+1j*rng.normal(size=8)
        perturb-=ref*np.vdot(ref,perturb);perturb/=np.linalg.norm(perturb)
        scale=rng.uniform(0.005,0.08)
        vec=np.sqrt(1-scale*scale)*ref+scale*perturb
        rho=0.999*np.outer(vec,vec.conj())+0.001*np.eye(8)/8
    refs=[];budget=0
    for site in range(3):
        order=[site]+[k for k in range(3) if k!=site]
        tensor=rho.reshape([2]*6).transpose(order+[k+3 for k in order]).reshape(2,4,2,4)
        red=np.einsum('atbt->ab',tensor)
        b=np.array([np.trace(red@m).real for m in [XN,YN,ZN]])
        delta=local_gap(b)
        budget+=1-np.trace(red@red).real+2*np.log(2)*delta
        refs.append(np.array(plus if b[1]>=0 else minus,dtype=complex).reshape(2))
    product=np.kron(np.kron(*refs[:2]),refs[2])
    fidelity=np.vdot(product,rho@product).real
    le(max(0,1-budget),fidelity,'network product fidelity budget')
    nontrivial_budgets += int(budget < 1)

result={'status':'passed','exact_checks':exact,'numerical_checks':numerical,'negative_controls':negative,
        'logical_cell_states':400,'three_node_states':200,'nontrivial_network_lower_bounds':nontrivial_budgets,'choi_erasure_minimum_eigenvalue':'-1/4',
        'maximum_isotropic_retention_at_zero_chirality':'1/2','oriented_volume_eigenvalues':['-sqrt(3)/4','sqrt(3)/4'],
        'python':platform.python_version(),'sympy':s.__version__,'numpy':np.__version__,
        'scope':'finite symbolic identities and seeded numerical inequalities; no Lean or universal-proof certification'}
print(json.dumps(result,indent=2))
