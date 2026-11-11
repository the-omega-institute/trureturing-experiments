#!/usr/bin/env python3
"""HS4 inter-cell exchange: exact identities and bounded spectral tests.

Uses only SymPy/NumPy/SciPy locally. No Lean, QPU or infinite-volume claim.
Every numerical assertion uses an explicit tolerance. Counts are assertions,
not individual entries of compared matrices. Run this file with no arguments.
"""
from __future__ import annotations
import itertools
import json
import platform
import numpy as np
import sympy as s
import scipy
from scipy import sparse

COUNTS={'exact_checks':0,'numerical_checks':0,'negative_controls':0}
def eq(a,b,label):
    d=a-b
    ok=all(s.simplify(v)==0 for v in d) if isinstance(d,s.MatrixBase) else s.simplify(d)==0
    if not ok: raise RuntimeError(label)
    COUNTS['exact_checks']+=1

def near(a,b,label,tol=2e-10):
    if np.max(np.abs(np.asarray(a)-np.asarray(b)))>tol: raise RuntimeError(label)
    COUNTS['numerical_checks']+=1

def le(a,b,label,tol=2e-10):
    if a>b+tol: raise RuntimeError(label)
    COUNTS['numerical_checks']+=1

def swap(n,i,j):
    entries={}
    for a in range(2**n):
        b=a
        if ((a>>(n-1-i))^(a>>(n-1-j)))&1:
            b^=(1<<(n-1-i))|(1<<(n-1-j))
        entries[b,a]=s.S.One
    return s.SparseMatrix(2**n,2**n,entries)

def kron(*args): return s.SparseMatrix(s.kronecker_product(*args))
def clean(m): return m.applyfunc(s.simplify)
I=s.eye(2); X=s.Matrix([[0,1],[1,0]]); Y=s.Matrix([[0,-s.I],[s.I,0]]); Z=s.diag(1,-1)
paulis=[I,X,Y,Z]
e0=s.zeros(16,1);e1=s.zeros(16,1)
for i in (3,12): e0[i]=1/s.sqrt(3)
for i in (5,6,9,10): e0[i]=-1/(2*s.sqrt(3))
for i,v in ((5,1),(6,-1),(9,-1),(10,1)): e1[i]=s.Rational(v,2)
E=s.SparseMatrix(e0.row_join(e1)); K=kron(E,E)
L=sum((swap(4,i,j) for i,j in itertools.combinations(range(4),2)),s.SparseMatrix(16,16,{}))
H0=(kron(L,s.eye(16))+kron(s.eye(16),L))/2
J={(i,a):kron(*[paulis[a+1]/2 if q==i else I for q in range(4)]) for i in range(4) for a in range(3)}
A={(i,k):3*I/4 if i==k else clean((2*E.T*swap(4,i,k)*E-I)/4) for i in range(4) for k in range(4)}
B=s.Matrix([[1,1,1],[1,-1,-1],[-1,1,-1],[-1,-1,1]])/2
C=s.eye(4)-s.ones(4)/4

eq(E.T*E,I,'code isometry');eq(L*E,s.zeros(16,2),'singlet closure')
eq(B.T*B,s.eye(3),'contrast isometry');eq(B*B.T,C,'contrast projection')
eq(L*(L-2*s.eye(16))*(L-6*s.eye(16)),s.zeros(16),'cell spectrum polynomial')
for i,a in itertools.product(range(4),range(3)):
    eq(E.T*J[i,a]*E,s.zeros(2),'one-spin compression')
    eq(L*J[i,a]*E,2*J[i,a]*E,'one-spin triplet image')
for i,k,a,b in itertools.product(range(4),range(4),range(3),range(3)):
    eq(E.T*J[i,a]*J[k,b]*E,(A[i,k]/3 if a==b else s.zeros(2)),'all second spin moments')
# Universal coefficient of logical Y in every ordered three-spin word.
left={(i,a):E.T*J[i,a] for i,a in J}
right={(k,b,m,c):J[k,b]*J[m,c]*E for k,b,m,c in itertools.product(range(4),range(3),range(4),range(3))}
for i,k,m in itertools.product(range(4),repeat=3):
    h=s.LeviCivita(i,k,m,next(iter(set(range(4))-{i,k,m}))) if len({i,k,m})==3 else 0
    for a,b,c in itertools.product(range(3),repeat=3):
        coefficient=s.trace(Y*left[i,a]*right[k,b,m,c])/2
        target=-s.sqrt(3)*h*s.LeviCivita(a,b,c)/24
        eq(coefficient,target,'all third-word chirality coefficients')
# General polynomial identity for an arbitrary, not necessarily symmetric, real T.
t=s.symbols('t0:16',real=True); T=s.Matrix(4,4,t)
cof=sum((-1)**(i+j)*T.minor_submatrix(i,j).det() for i in range(4) for j in range(4))
eq(s.expand(cof-4*(B.T*T*B).det()),0,'all-variable determinant/cofactor identity')
I256=s.SparseMatrix(s.eye(256))
links={(i,j):swap(8,i,4+j)/2-I256/4 for i,j in itertools.product(range(4),repeat=2)}
def make_v(T):
    return sum((T[i,j]*links[i,j] for i,j in itertools.product(range(4),repeat=2) if T[i,j]),s.SparseMatrix(256,256,{}))
matrices=[s.eye(4),s.diag(1,1,1,0),s.diag(1,1,0,0),s.diag(1,0,0,0),s.ones(4),
          s.Matrix([[0,1,0,0],[1,0,0,0],[0,0,1,0],[0,0,0,1]]),
          s.Matrix([[1,2,0,-1],[0,1,3,1],[2,-1,1,0],[0,2,1,1]]),
          B*s.diag(2,-1,0)*B.T,
          s.ones(4,1)*s.Matrix([[1,2,3,4]])+s.Matrix([2,-1,1,3])*s.ones(1,4)]
for T in matrices:
    V=make_v(T);W=V*K;M2=clean(W.T*W);M3=clean(W.T*V*W)
    eq(K.T*W,s.zeros(4),'first order zero')
    eq(H0*W,2*W,'two-cell fixed denominator')
    eq(s.trace(M2)/4,s.trace((B.T*T*B).T*(B.T*T*B))/3,'mean virtual shift')
    predicted=sum((T[i,j]*T[k,l]*kron(A[i,k],A[j,l])/3 for i,j,k,l in itertools.product(range(4),repeat=4) if T[i,j]*T[k,l]),s.zeros(4))
    eq(M2,predicted,'full second order contraction')
    eq(s.trace(kron(Y,Y)*M2)/4,0,'second order YY absence')
    eq(s.trace(kron(Y,Y)*M3)/16,3*(B.T*T*B).det()/16,'third order determinant coefficient')
    if B.T*T*B==s.zeros(3): eq(W,s.zeros(256,4),'zero contrast exact darkness')
V=make_v(s.eye(4));W=V*K
M2=clean(W.T*W);M3=clean(W.T*V*W)
eq(-M2/2,-kron(I,I)/2-(kron(X,X)+kron(Z,Z))/4,'matched second order')
eq(M3/4,-3*(kron(I,I)+kron(X,X)+kron(Z,Z)-kron(Y,Y))/16,'matched third order')
phi=s.Matrix([1,0,0,1])/s.sqrt(2); pp=phi*phi.T
for i,j in itertools.combinations(range(4),2):
    S=clean(E.T*swap(4,i,j)*E)
    eq(kron(S,S)*phi,phi,'paired permutation invariant')
# The order-three Fourier measurement, with exact algebraic roots of unity.
w=(-1+s.I*s.sqrt(3))/2
cyc=swap(4,0,1)*swap(4,1,2)
eq(cyc**3,s.eye(16),'physical permutation order')
eq(E.T*cyc*E,-I/2+s.I*s.sqrt(3)*Y/2,'three-cycle logical representation')
for k,target in [(0,s.zeros(2)),(1,(I+Y)/2),(2,(I-Y)/2)]:
    F=sum((w**(-k*l)*cyc**l for l in range(3)),s.zeros(16))/3
    eq(E.T*F*E,target,'character projectors')

# Bounded numerical tests are separate from the exact assertions above.
HN=np.array(H0,dtype=float);VN=np.array(V,dtype=float);KN=np.array(K,dtype=float)
M2N=np.array(M2,dtype=float);M3N=np.array(M3,dtype=float)
base_evals,base_vecs=np.linalg.eigh(HN)
Qbasis=base_vecs[:,base_evals>0.5]
qn=Qbasis.T@HN@Qbasis
pn=KN@np.array(phi,dtype=float).reshape(4)
records=[]
for g in [0.001,0.002,0.005,0.01,0.02,0.04,0.08]:
    U=g*VN; evals,vecs=np.linalg.eigh(HN+U); vnorm=3*g
    eff=-g*g*M2N/2+g**3*M3N/4
    mu=np.linalg.eigvalsh(eff);delta=vnorm**4
    le(np.max(np.abs(evals[:4]-mu)),delta,'four eigenvalues with uniform remainder')
    for lam in evals[:4]:
        WN=Qbasis.T@U@KN
        FN=-WN.T@np.linalg.solve(qn+Qbasis.T@U@Qbasis-lam*np.eye(252),WN)
        le(np.linalg.norm(FN-eff,2),delta,'on-shell Schur complement bound')
    lower=g*g/2+3*g**3/4-2*delta
    if lower>0:
        le(lower,evals[1]-evals[0],'certified full-system gap')
        fidelity=float(abs(pn@vecs[:,0])**2)
        le(1/(1+(vnorm/(1-vnorm))**2),fidelity,'certified encoded Bell fidelity')
        projected=KN.T@vecs[:,0];projected/=np.linalg.norm(projected)
        near(np.outer(projected,projected),np.array(pp,dtype=float),'ground projected exactly on Bell line')
    records.append({'g_over_Delta':g,'energy0':float(evals[0]),'gap':float(evals[1]-evals[0]),
                    'maximum_third_order_energy_error':float(np.max(np.abs(evals[:4]-mu))),
                    'proved_error_bound':float(delta),'encoded_Bell_overlap':float(abs(pn@vecs[:,0])**2)})
rng=np.random.default_rng(614203)
for _ in range(8):
    T=s.Matrix(rng.integers(-2,3,size=(4,4)).tolist());U0=make_v(T)
    U=np.array(U0,dtype=float);U*=0.08/max(1,np.linalg.norm(U,2))
    vnorm=float(np.linalg.norm(U,2));WN=U@KN
    eff=-WN.T@WN/2+WN.T@U@WN/4
    evals=np.linalg.eigvalsh(HN+U)
    le(float(np.max(np.abs(evals[:4]-np.linalg.eigvalsh(eff)))),vnorm**4,'random coupling spectral bound')

# A three-cell path: no triangle term in the cubic coefficient.
def numeric_swap(n,i,j):
    a=np.arange(2**n);b=a.copy(); mask=(((a>>(n-1-i))^(a>>(n-1-j)))&1).astype(bool)
    b[mask]^=(1<<(n-1-i))|(1<<(n-1-j))
    return sparse.csr_matrix((np.ones(2**n),(b,a)),shape=(2**n,2**n))
K3=np.kron(np.kron(np.array(E,dtype=float),np.array(E,dtype=float)),np.array(E,dtype=float))
I4096=sparse.eye(4096,format='csr');VT=sparse.csr_matrix((4096,4096))
for a,b in [(0,1),(1,2)]:
    for i in range(4): VT+=numeric_swap(12,4*a+i,4*b+i)/2-I4096/4
W3=VT@K3
near(-W3.T@W3/2,np.kron(-M2N/2,np.eye(2))+np.kron(np.eye(2),-M2N/2),'tree second coefficient')
near(W3.T@(VT@W3)/4,np.kron(M3N/4,np.eye(2))+np.kron(np.eye(2),M3N/4),'tree third coefficient')
# Three-cell cycles have an additional linked term, computed independently.
triangle_exact=sum((kron(A[i,k],A[i,j],A[j,k])/6 for i,j,k in itertools.product(range(4),repeat=3)),s.zeros(8))
triangle_closed=kron(I,I,I)/2+(kron(X,X,I)+kron(X,I,X)+kron(I,X,X)+kron(Z,Z,I)+kron(Z,I,Z)+kron(I,Z,Z))/4+(kron(Z,Z,Z)-kron(X,X,Z)-kron(X,Z,X)-kron(Z,X,X))/8
eq(triangle_exact,triangle_closed,'full triangle linked coefficient')
edge_ops=[]
for a,b in [(0,1),(1,2),(0,2)]:
    U=sparse.csr_matrix((4096,4096))
    for i in range(4): U+=numeric_swap(12,4*a+i,4*b+i)/2-I4096/4
    edge_ops.append(U)
for weights in [(1,1,1),(2,-1,3)]:
    allv=sum((g*U for g,U in zip(weights,edge_ops)),sparse.csr_matrix((4096,4096)))
    allw=allv@K3
    individual=sum((g**3*(U@K3).T@(U@(U@K3))/4 for g,U in zip(weights,edge_ops)),np.zeros((8,8)))
    near(allw.T@(allv@allw)/4,individual+np.prod(weights)*np.array(triangle_closed,dtype=float),'triangle correction on a weighted cycle')
# Concrete false claims must be rejected, not admitted by assert removal.
for wrong in [s.Rational(3,16)==0, s.det(B.T*s.diag(1,1,0,0)*B)!=0,
              clean(E.T*cyc*E)==I]:
    if bool(wrong): raise RuntimeError('negative control accepted')
    COUNTS['negative_controls']+=1
print(json.dumps({'status':'passed',**COUNTS,'exact_coupling_matrices':len(matrices),
    'random_spectral_matrices':8,'matched_spectral_cases':records,
    'third_order_YY_coefficient':'3 det(B^T T B)/(16 Delta^2)',
    'at_g_over_Delta_0_01':{'proved_gap_lower_bound':'4913/100000000 * Delta',
                           'proved_encoded_Bell_fidelity_lower_bound':'9409/9418'},
    'python':platform.python_version(),'numpy':np.__version__,'sympy':s.__version__,'scipy':scipy.__version__,
    'scope':'exact finite identities and bounded numerical tests; no Lean, hardware, or thermodynamic-limit certificate'},indent=2))
