#!/usr/bin/env python3
"""One extra phase address: native routing, full-circuit compilation and witness readout.
Finite exact and numerical evidence for the specified ideal no11 model.
"""
from __future__ import annotations
import collections
import hashlib
import itertools
import json
import math
import platform
from pathlib import Path
import numpy as np
import sympy as sp

COUNTS = collections.Counter()
I = sp.eye(2)
X = sp.Matrix([[0, 1], [1, 0]])
Y = sp.Matrix([[0, -sp.I], [sp.I, 0]])
Z = sp.diag(1, -1)
P0, P1 = sp.diag(1, 0), sp.diag(0, 1)

def eq(a, b, label):
    d = a-b
    entries = list(d) if isinstance(d, sp.MatrixBase) else [d]
    if not all(sp.simplify(x) == 0 for x in entries):
        raise RuntimeError(label)
    COUNTS['symbolic_equalities'] += 1

def check(ok, label, category='structural_assertions'):
    if not bool(ok):
        raise RuntimeError(label)
    COUNTS[category] += 1

def near(a, b, label, tol=2e-10):
    check(np.max(np.abs(np.asarray(a)-np.asarray(b))) <= tol,
          label, 'floating_comparisons')

def words(n):
    return [w for w in itertools.product((0, 1), repeat=n)
            if all(w[k]*w[k+1] == 0 for k in range(n-1))]

def edges(n):
    ws = words(n); ix = {w:i for i,w in enumerate(ws)}
    ans = []
    for k in range(n):
        out = []
        for a,w in enumerate(ws):
            if w[k] == 0 and (k == 0 or not w[k-1]) and (k == n-1 or not w[k+1]):
                v = list(w); v[k] = 1; out.append((a,ix[tuple(v)]))
        ans.append(out)
    return ans

def encoding(m):
    ws=words(2*m-1); ix={w:i for i,w in enumerate(ws)}
    E=sp.zeros(len(ws),2**m)
    for a,bits in enumerate(itertools.product((0,1),repeat=m)):
        w=tuple(v for j,b in enumerate(bits) for v in ([b,0] if j<m-1 else [b]))
        E[ix[w],a]=1
    return E

def one(m,q,gate):
    return sp.kronecker_product(*[gate if j==q else I for j in range(m)])

def perm_gate(m,c,t,kind):
    ans=sp.zeros(2**m)
    for a in range(2**m):
        bits=[(a>>(m-1-j))&1 for j in range(m)]
        amp=1; out=bits.copy()
        if kind=='SWAP':out[c],out[t]=out[t],out[c]
        elif kind=='CX' and bits[c]:out[t]^=1
        elif kind=='CY' and bits[c]:
            amp=sp.I if bits[t]==0 else -sp.I;out[t]^=1
        elif kind=='CZ' and bits[c] and bits[t]:amp=-1
        b=sum(v<<(m-1-j) for j,v in enumerate(out));ans[b,a]=amp
    return ans

class Compiler:
    """Pulse list only. All beta=pi/2 pulses use the selected data address."""
    def __init__(self,m,reference):
        if m<1 or not 0<=reference<m:raise ValueError('invalid register/reference')
        self.m=m;self.reference=reference;self.pulses=[];self.phase=1+0j
    def pulse(self,site,beta,t):
        if beta not in (0,sp.pi/2):raise ValueError('unsupported quadrature')
        if beta and site!=2*self.reference:raise ValueError('second phase outside reference')
        self.pulses.append((site,beta,t))
    def rx(self,q,t):self.pulse(2*q,0,t)
    def cz_adj(self,a,b):
        if abs(a-b)!=1:raise ValueError('nonadjacent CZ')
        q=min(a,b)
        self.rx(q,sp.pi/2);self.rx(q+1,sp.pi/2)
        self.pulse(2*q+1,0,sp.pi)
        self.rx(q+1,-sp.pi/2);self.rx(q,-sp.pi/2)
    def cy_adj(self,c,t):
        self.rx(t,sp.pi/4);self.cz_adj(c,t);self.rx(t,-sp.pi/4)
    def swap_adj(self,a,b):
        if abs(a-b)!=1:raise ValueError('nonadjacent SWAP')
        self.cy_adj(a,b);self.cy_adj(b,a);self.cy_adj(a,b)
    def route(self,a,b):
        steps=[];p=a
        while p!=b:
            q=p+(1 if b>p else -1)
            self.swap_adj(p,q);steps.append((p,q));p=q
        return steps
    def ry(self,q,t):
        steps=self.route(q,self.reference)
        self.pulse(2*self.reference,sp.pi/2,t)
        for a,b in reversed(steps):self.swap_adj(a,b)
    def cz(self,a,b):
        if a==b:raise ValueError('distinct operands required')
        p=a;step=1 if b>a else -1;steps=[]
        while abs(p-b)>1:
            self.swap_adj(p,p+step);steps.append((p,p+step));p+=step
        self.cz_adj(p,b)
        for c,t in reversed(steps):self.swap_adj(c,t)
    def h_exact_phase(self,q,inverse=False):
        if inverse:
            self.rx(q,-sp.pi/2);self.ry(q,-sp.pi/4)
        else:
            self.ry(q,sp.pi/4);self.rx(q,sp.pi/2)
    def h(self,q):
        self.h_exact_phase(q);self.phase*=-1j
    def phase_gate(self,q,t):
        self.rx(q,-sp.pi/4);self.ry(q,t/2);self.rx(q,sp.pi/4)
        self.phase*=np.exp(-.5j*float(t))
    def cx(self,c,t):
        self.h_exact_phase(t,True);self.cz(c,t);self.h_exact_phase(t)
    def execute(self,ops):
        for op in ops:
            if op[0]=='H':self.h(op[1])
            elif op[0]=='P':self.phase_gate(op[1],op[2])
            elif op[0]=='CX':self.cx(op[1],op[2])
            else:raise ValueError(op)

SP_CACHE={}
def symbolic_matrix(n,pulses):
    es=edges(n);d=len(words(n));U=sp.eye(d)
    for site,beta,t in pulses:
        key=(n,site,beta,t)
        if key not in SP_CACHE:
            g=sp.eye(d)
            for a,b in es[site]:
                g[a,a]=g[b,b]=sp.cos(t)
                g[a,b]=-sp.I*sp.sin(t)*sp.exp(-sp.I*beta)
                g[b,a]=-sp.I*sp.sin(t)*sp.exp(sp.I*beta)
            SP_CACHE[key]=g
        U=(SP_CACHE[key]*U).applyfunc(sp.expand)
    return U

def numerical_action(n,pulses,state):
    es=edges(n);out=np.asarray(state,dtype=complex).copy()
    for site,beta,t in pulses:
        beta=float(beta);t=float(t);c=np.cos(t);s=np.sin(t)
        aa=np.array([a for a,b in es[site]],dtype=int)
        bb=np.array([b for a,b in es[site]],dtype=int)
        a=out[aa].copy();b=out[bb].copy()
        out[aa]=c*a-1j*s*np.exp(-1j*beta)*b
        out[bb]=c*b-1j*s*np.exp(1j*beta)*a
    return out

# Exact full-space compressed actions of the zero-phase routing primitives.
for m in [2,3]:
    E=encoding(m);n=2*m-1
    for a in range(m-1):
        for c,t in [(a,a+1),(a+1,a)]:
            comp=Compiler(m,0);comp.cy_adj(c,t)
            check(len(comp.pulses)==7 and all(beta==0 for _,beta,_ in comp.pulses),'CY budget/menu')
            eq(symbolic_matrix(n,comp.pulses)*E,E*perm_gate(m,c,t,'CY'),'native CY')
        comp=Compiler(m,0);comp.swap_adj(a,a+1)
        check(len(comp.pulses)==21 and all(beta==0 for _,beta,_ in comp.pulses),'SWAP budget/menu')
        eq(symbolic_matrix(n,comp.pulses)*E,E*perm_gate(m,a,a+1,'SWAP'),'native SWAP')

# Verify arbitrary-angle reference transport on the full five-dimensional cell.
theta=sp.symbols('theta',real=True)
for h in [0,1]:
    for q in [0,1]:
        comp=Compiler(2,h);comp.ry(q,theta)
        eq(symbolic_matrix(3,comp.pulses)*encoding(2),
           encoding(2)*one(2,q,sp.cos(theta)*I-sp.I*sp.sin(theta)*Y),'symbolic remote Y')
        check(len(comp.pulses)==42*abs(q-h)+1,'remote Y pulse budget')

# Latest dev C-121 PVM on ALL five modes, including coherences with the apex.
E=encoding(2);ws=words(3);d=len(ws);apex=ws.index((0,1,0))
Pb=E*E.T;ap=sp.eye(d)-Pb
Ks=[E*sp.kronecker_product(X,Z)*E.T,E*sp.kronecker_product(Z,X)*E.T]
measurements=[]
for setting in [0,1]:
    comp=Compiler(2,0)
    if setting==1:comp.swap_adj(0,1)
    comp.cz_adj(0,1);comp.h_exact_phase(0,True)
    U=symbolic_matrix(3,comp.pulses)
    check(len(comp.pulses)==[7,28][setting],'witness setting pulse count')
    check(sum(beta!=0 for _,beta,_ in comp.pulses)==1,'one quadrature call per witness')
    for outcome in [1,-1,0]:
        F=sp.zeros(d)
        for j,w in enumerate(ws):
            if (outcome==0 and j==apex) or (outcome!=0 and j!=apex and (-1)**w[0]==outcome):F[j,j]=1
        target=ap if outcome==0 else (Pb+outcome*Ks[setting])/2
        eq(U.H*F*U,target,'complete five-mode PVM effect')
    measurements.append((comp,U))

# Exact source example and sharp range-four TV drift bound example.
e0=E*sp.Matrix([1,1,0,0])/sp.sqrt(2)
e1=E*sp.Matrix([0,0,1,-1])/sp.sqrt(2)
gamma,r=sp.symbols('gamma r',real=True)
rho=(1-r)*ap+r*(e0*e0.H+e1*e1.H+gamma*(e0*e1.H+e1*e0.H))/2
eq(sp.trace(rho*(Ks[0]+Ks[1]-Pb)),r*gamma,'source score')
eps=sp.symbols('eps',real=True)
p=sp.Matrix([sp.Rational(3,4),sp.Rational(1,4),0]);q=p+sp.Matrix([eps,-eps,0])
score=sp.Matrix([1,-3,0])
eq((score.T*p)[0],0,'boundary null score')
eq((score.T*q)[0],4*eps,'TV error gives exact range-four drift')
eq((score.T*q)[0]-4*eps,0,'corrected drift')

# Numerical all-input checks for one reference at every data address.
rng=np.random.default_rng(20261016)
max_route_residual=0.
for m in [1,2,3,4]:
    En=np.array(encoding(m),complex);n=2*m-1
    for h in range(m):
        for target in range(m):
            t=float(rng.uniform(-1,1));comp=Compiler(m,h);comp.ry(target,t)
            output=numerical_action(n,comp.pulses,En)
            wanted=En@np.array(one(m,target,np.cos(t)*I-1j*np.sin(t)*Y),complex)
            near(output,wanted,'remote single-reference action')
            max_route_residual=max(max_route_residual,float(np.linalg.norm(output-wanted,2)))
            check(all(beta==0 or site==2*h for site,beta,_ in comp.pulses),'reference location enforced')

# A complete known five-qubit order-probe, with the original small repeated powers.
def ccx(a,b,t):
    p=sp.pi/4
    return [('H',t),('CX',b,t),('P',t,-p),('CX',a,t),('P',t,p),('CX',b,t),
            ('P',t,-p),('CX',a,t),('P',b,p),('P',t,p),('H',t),('CX',a,b),('P',a,p),('P',b,-p),('CX',a,b)]
def cp(a,b,t):
    return [('P',a,t/2),('P',b,t/2),('CX',a,b),('P',b,-t/2),('CX',a,b)]
def sw(a,b):return [('CX',a,b),('CX',b,a),('CX',a,b)]
ops=[('H',i) for i in range(3)]
for clock in range(3):
    for _ in range(2**(2-clock)):
        ops+=ccx(clock,4,3)+ccx(clock,3,4)
forward=[]
for a in range(3):
    forward.append(('H',a))
    for b in range(a+1,3):forward+=cp(b,a,sp.pi/2**(b-a))
forward+=sw(0,2)
ops += [('P',op[1],-op[2]) if op[0]=='P' else op for op in reversed(forward)]

m=5;h=2;comp=Compiler(m,h);comp.execute(ops)
En=np.array(encoding(m),complex);native=numerical_action(9,comp.pulses,En)
# Independent matrix: uniform clock preparation, controlled Fibonacci powers, inverse QFT.
had=np.array((X+Z)/sp.sqrt(2),complex);prep=np.kron(np.kron(had,had),had)
U=np.zeros((4,4),complex)
for x in range(2):
    for y in range(2):U[2*(x^y)+x,2*x+y]=1
controlled=np.zeros((32,32),complex)
for a in range(8):controlled[4*a:4*a+4,4*a:4*a+4]=np.linalg.matrix_power(U,a)
F=np.exp(-2j*np.pi*np.outer(np.arange(8),np.arange(8))/8)/np.sqrt(8)
ideal=np.kron(F,np.eye(4))@controlled@np.kron(prep,np.eye(4))
near(native,comp.phase*En@ideal,'entire single-reference order-probe')
res=float(np.linalg.norm(native-comp.phase*En@ideal,2))
leak=float(np.linalg.norm(native-En@(En.conj().T@native),2))
check(all(beta==0 or site==2*h for site,beta,_ in comp.pulses),'full circuit one phase address')

# Perturbed measurement comparisons: unitary error controls the actual output law.
max_tv=0.
for setting,(co,Uexact) in enumerate(measurements):
    Un=np.array(Uexact,complex)
    for _ in range(12):
        z=rng.normal(size=(5,5))+1j*rng.normal(size=(5,5));state=z@z.conj().T;state/=np.trace(state)
        a=rng.normal(size=(5,5))+1j*rng.normal(size=(5,5));a=(a+a.conj().T)/2
        vals,vecs=np.linalg.eigh(a);a/=max(abs(vals));vals/=max(abs(vals))
        noise=(vecs*np.exp(-1j*.002*vals))@vecs.conj().T
        noisy=noise@Un;err=np.linalg.norm(noisy-Un,2)
        p=np.real(np.diag(Un@state@Un.conj().T));q=np.real(np.diag(noisy@state@noisy.conj().T))
        tv=float(np.sum(abs(p-q))/2);max_tv=max(max_tv,tv)
        check(tv<=err+1e-12,'unitary error versus word total variation','floating_comparisons')
# Concrete rejected claims. These are sanity controls, not independent proof review.
CY=perm_gate(2,0,1,'CY');CX=perm_gate(2,0,1,'CX')
wrong=[CY==CX, 4*sp.Rational(1,100)<=0,
       all(beta==0 for _,beta,_ in comp.pulses)]
for statement in wrong:
    check(not bool(statement),'negative control unexpectedly accepted','negative_controls')

serialized=[[site+1,float(beta),float(t)] for site,beta,t in comp.pulses]
pulsebytes=json.dumps(serialized,separators=(',',':')).encode()
Path(__file__).with_name('single_reference_pulses.json').write_bytes(pulsebytes+b'\n')
result={'status':'passed',**dict(COUNTS),
    'source_gates':len(ops),'native_pulses':len(comp.pulses),
    'logical_qubits':5,'native_sites':9,'reference_site_one_based':2*h+1,
    'second_quadrature_pulses':sum(beta!=0 for _,beta,_ in comp.pulses),
    'second_quadrature_sites_one_based':sorted({site+1 for site,beta,_ in comp.pulses if beta!=0}),
    'full_encoded_operator_residual':res,'final_separator_leakage_norm':leak,
    'max_remote_rotation_residual':max_route_residual,
    'witness_pulses_by_setting':[len(x[0].pulses) for x in measurements],
    'max_sampled_measurement_tv':max_tv,
    'pulse_list_sha256':hashlib.sha256(pulsebytes+b'\n').hexdigest(),
    'python':platform.python_version(),'sympy':sp.__version__,'numpy':np.__version__,
    'scope':'ideal native-model symbolic identities and finite floating checks; no Lean, QPU, optimal-pulse or quantum-speedup certificate'}
print(json.dumps(result,indent=2))
