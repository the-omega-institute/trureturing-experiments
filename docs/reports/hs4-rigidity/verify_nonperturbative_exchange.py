#!/usr/bin/env python3
"""Exact identities and bounded tests for nonperturbative HS4 cell exchange.

The general positivity, uniqueness, majorization and graph-bound proofs are
in the theory. Numerical checks are not universal or interval certificates.
"""
from __future__ import annotations
import itertools
import json
import platform
import numpy as np
import scipy
from scipy import sparse
from scipy.sparse.linalg import eigsh
import sympy as s

counts = dict(exact_checks=0, numerical_checks=0, negative_controls=0)

def eq(a, b, name):
    v = a-b
    good = all(s.simplify(x)==0 for x in v) if isinstance(v,s.MatrixBase) else s.simplify(v)==0
    if not good: raise RuntimeError(name)
    counts['exact_checks'] += 1

def near(a,b,name,tol=2e-10):
    if np.max(np.abs(np.asarray(a)-np.asarray(b)))>tol: raise RuntimeError(name)
    counts['numerical_checks'] += 1

def le(a,b,name,tol=2e-10):
    if float(a)>float(b)+tol: raise RuntimeError(name)
    counts['numerical_checks'] += 1

I=s.eye(2); X=s.Matrix([[0,1],[1,0]]); Y=s.Matrix([[0,-s.I],[s.I,0]]); Z=s.diag(1,-1)
pauli=[X,Y,Z]

def kron(*xs):return s.kronecker_product(*xs)

def swap(n,i,j):
    d=s.zeros(2**n)
    for a in range(2**n):
        b=a
        if ((a>>(n-1-i))^(a>>(n-1-j)))&1:b^=(1<<(n-1-i))|(1<<(n-1-j))
        d[b,a]=1
    return d

L=sum((swap(4,i,j) for i,j in itertools.combinations(range(4),2)),s.zeros(16))
I16=s.eye(16)
P=[(L-2*I16)*(L-6*I16)/12,L*(6*I16-L)/8,L*(L-2*I16)/24]
dims=[2,9,5]
J=[kron(*[p/2 if k==i else I for k in range(4)]) for i in range(4) for p in pauli]
R=kron(Y,Y,Y,Y)

eq(R.T*R,I16,'time-reversal linear part unitary')
eq(R.T,R,'four-spin reversal symmetric')
eq(R*L.T*R,L,'cell time-reversal symmetry')
for A in J:eq(R*A.T*R,-A,'all local spin reversal laws')
for i in range(3):
    eq(P[i].T,P[i],'projector real symmetric')
    eq(s.trace(P[i]),dims[i],'projector rank')
    eq(L*P[i],[0,2,6][i]*P[i],'projector spin value')
    for j in range(3):eq(P[i]*P[j],P[i] if i==j else s.zeros(16),'projector orthogonality')
eq(sum(P,s.zeros(16)),I16,'complete spin projectors')

def Q(A):return sum((j*A*j for j in J),s.zeros(16))
eq(Q(I16),3*I16,'sandwich identity')
eq(Q(L),L+6*I16,'sandwich first moment')
eq(Q(L*L),-L*L+14*L+12*I16,'sandwich second moment')
Qcols=[2*P[1]/3,3*P[0]+3*P[1]/2+3*P[2]/2,5*P[1]/6+3*P[2]/2]
for i in range(3):eq(Q(P[i]),Qcols[i],'all-sector sandwich closure')
D,g,z=s.symbols('Delta g z',real=True)
M=s.Matrix([[0,-s.sqrt(2)*g,0],[-s.sqrt(2)*g,2*D-3*g/2,-s.sqrt(5)*g/2],[0,-s.sqrt(5)*g/2,6*D-3*g/2]])
Mcalc=s.Matrix(3,3,lambda i,j:s.simplify(s.trace(P[i]*(D*(L*P[j]+P[j]*L)/2-g*Q(P[j])))/s.sqrt(dims[i]*dims[j])))
eq(Mcalc,M,'exact three-by-three block')
poly=z**3+(3*g-8*D)*z**2+(12*D**2-12*D*g-g**2)*z+12*D*g**2-3*g**3
eq((z*s.eye(3)-M).det(),poly,'ground energy cubic')
strong=s.Matrix([s.sqrt(2),3,s.sqrt(5)])/4
eq(M.subs(D,0)*strong,-3*g*strong,'strong-coupling endpoint')
e0=s.zeros(16,1);e1=s.zeros(16,1)
for i in [3,12]:e0[i]=1/s.sqrt(3)
for i in [5,6,9,10]:e0[i]=-1/(2*s.sqrt(3))
for i,v in [(5,1),(6,-1),(9,-1),(10,1)]:e1[i]=s.Rational(v,2)
E=e0.row_join(e1)
vec=lambda A:s.Matrix(list(A))
Xi=[vec(P[i]*R.T)/s.sqrt(dims[i]) for i in range(3)]
U=s.Matrix.hstack(*Xi)
eq(U.H*U,s.eye(3),'physical block embedding')
eq(R*E,E,'real singlets fixed by reversal')
Bell=(kron(e0,e0)+kron(e1,e1))/s.sqrt(2)
eq(Xi[0],Bell,'HS handed-pair Bell identity')
# Check the entire microscopic action on the invariant subspace, separately.
I256=s.eye(256)
H0=(kron(L,I16)+kron(I16,L))/2
V=sum((swap(8,i,4+i)/2-I256/4 for i in range(4)),s.zeros(256))
eq(H0*U,U*M.subs({D:1,g:0}),'full eight-spin free action')
eq(V*U,U*M.subs({D:0,g:1}),'full eight-spin exchange action')
# A parameter-independent one-way instrument in sector-adapted bases.
# Its Bob completion resets the unused basis states; selected inputs never use it.
Aeffects=s.zeros(16); branches=[];offset=0
for j,d in enumerate(dims):
    for k,l in itertools.combinations(range(offset,offset+d),2):
        A=s.zeros(2,16); A[0,k]=1/s.sqrt(d-1);A[1,l]=1/s.sqrt(d-1)
        B=s.zeros(2,16);B[0,k]=1;B[1,l]=1
        Aeffects+=A.H*A
        eq(B.H*B+s.diag(*[0 if q in (k,l) else 1 for q in range(16)]),I16,'Bob trace-preserving completion')
        branches.append((j,A,B))
    offset+=d
eq(Aeffects,I16,'47-outcome Alice normalization')
canonical=[];off=0
for d in dims:
    v=s.zeros(256,1)
    for k in range(off,off+d):v[16*k+k]=1/s.sqrt(d)
    canonical.append(v);off+=d
phi=s.Matrix([1,0,0,1])/s.sqrt(2)
for j in range(3):
    for k in range(3):
        out=s.zeros(4)
        for sector,A,B in branches:
            cj=s.diag(*[1/s.sqrt(dims[j]) if sum(dims[:j])<=q<sum(dims[:j+1]) else 0 for q in range(16)])
            ck=s.diag(*[1/s.sqrt(dims[k]) if sum(dims[:k])<=q<sum(dims[:k+1]) else 0 for q in range(16)])
            vj=vec(A*cj*B.T);vk=vec(A*ck*B.T)
            out+=vj*vk.H
        eq(out,phi*phi.H if j==k else s.zeros(4),'all nine sector-matrix decoder actions')

MN=s.lambdify((D,g),M,'numpy')
UN=np.array(U.evalf(),complex);H0N=np.array(H0,float);VN=np.array(V,float)
PN=[np.array(x,float) for x in P];RN=np.array(R,float)
reference=[];worst_energy=worst_vector=0.0
# Every floating comparison uses a scaled Hamiltonian with max(Delta,g)=1.
for ratio in np.geomspace(1e-3,1e3,37):
    scale=max(1.0,ratio);dd=1/scale;gg=ratio/scale
    vals,vs=np.linalg.eigh(MN(dd,gg));alpha=vs[:,0];alpha*=np.sign(alpha[0]);c=alpha/np.sqrt(dims)
    true,states=np.linalg.eigh(dd*H0N+gg*VN);psi=states[:,0];pred=UN@alpha
    energyerr=abs(true[0]-vals[0]);vectorerr=1-abs(np.vdot(pred,psi))**2
    worst_energy=max(worst_energy,energyerr);worst_vector=max(worst_vector,abs(vectorerr))
    near(true[0],vals[0],'full ground energy versus exact block')
    le(vectorerr,0,'full ground vector versus exact block')
    le(1e-12,true[1]-true[0],'sampled positive full gap')
    C=(psi.reshape(16,16)@RN).astype(complex);C*=np.exp(-1j*np.angle(np.trace(C)))
    near(C,C.conj().T,'ground coefficient Hermitian')
    near(C,sum(c[i]*PN[i] for i in range(3)),'positive ground coefficient reconstruction')
    lam=np.repeat(c*c,dims)
    near(np.linalg.eigvalsh(C@C.conj().T),np.sort(lam),'complete Schmidt spectrum')
    le(0,c[2],'sampled strict Schmidt positivity')
    le(c[2],c[1],'ordered smallest Schmidt value');le(c[1],c[0],'ordered largest Schmidt value')
    p0=alpha[0]**2;le(1/8,p0,'heralded probability lower bound');le(p0,1,'heralded probability upper bound')
    entropy=float(-sum(lam*np.log2(lam)));le(1,entropy,'one-ebit lower bound');le(entropy,4,'four-ebit upper bound')
    if reference:
        prior=reference[-1]['schmidt_spectrum']
        le(float(np.max(np.cumsum(lam)[:-1]-np.cumsum(prior)[:-1])),0,'sampled majorization chain')
    t=1/ratio;ev,av=np.linalg.eigh(MN(t,1));a=av[:,0];a*=np.sign(a[0]);ep=2*a[1]**2+6*a[2]**2
    le(0,ep,'Hellmann-Feynman derivative positive');le(ep,6,'Hellmann-Feynman derivative below six')
    near(c[1]/c[0],-ev[0]/3,'first Schmidt ratio formula')
    near(c[2]/c[1],1.5/(6*t-1.5-ev[0]),'second Schmidt ratio formula')
    reference.append(dict(g_over_Delta=float(ratio),energy_over_Delta=float(vals[0]*scale),
        full_gap_over_Delta=float((true[1]-true[0])*scale),entropy_bits=entropy,
        singlet_success=float(p0),schmidt_spectrum=lam.tolist()))
# Independently disordered positive matched links retain uniqueness/full rank.
rng=np.random.default_rng(2026101407)
linksN=[np.array(swap(8,i,4+i)/2-I256/4,float) for i in range(4)]
disordered_min_rank_value=1.0
for _ in range(12):
    weights=rng.uniform(.25,3.0,4)
    hh=H0N+sum(w*b for w,b in zip(weights,linksN));ev,vs=np.linalg.eigh(hh)
    cm=vs[:,0].reshape(16,16)@RN;cm*=np.sign(np.trace(cm))
    near(cm,cm.T,'unequal-link coefficient Hermitian')
    le(1e-10,np.linalg.eigvalsh(cm)[0],'unequal-link full rank')
    le(1e-10,ev[1]-ev[0],'unequal-link unique ground sample')
    disordered_min_rank_value=min(disordered_min_rank_value,float(np.linalg.eigvalsh(cm)[0]))
# Three-cell complete Hamiltonians for graph energy brackets.
def sp_swap(n,i,j):
    a=np.arange(2**n);b=a.copy();mask=(((a>>(n-1-i))^(a>>(n-1-j)))&1).astype(bool)
    b[mask]^=(1<<(n-1-i))|(1<<(n-1-j))
    return sparse.csr_matrix((np.ones(2**n),(b,a)),shape=(2**n,2**n))

def graph_h(n,edges,gg):
    eye=sparse.eye(2**(4*n),format='csr');hh=sparse.csr_matrix(eye.shape)
    for a in range(n):
        for i,j in itertools.combinations(range(4),2):hh+=sp_swap(4*n,4*a+i,4*a+j)/2
    for a,b in edges:
        for i in range(4):hh+=gg*(sp_swap(4*n,4*a+i,4*b+i)/2-eye/4)
    return hh
brackets=[]
for name,edges in [('path',[(0,1),(1,2)]),('triangle',[(0,1),(1,2),(0,2)])]:
    for gg in [.3,2.0]:
        hh=graph_h(3,edges,gg)
        # All terms conserve total z magnetization; a total-spin eigenstate has an m=0 representative.
        idx=np.array([j for j in range(4096) if int(j).bit_count()==6]);hz=hh[idx,:][:,idx]
        vv,ss=eigsh(hz,k=1,which='SA',tol=1e-12,maxiter=20000,v0=rng.normal(size=len(idx)))
        value=float(vv[0]);near(hz@ss[:,0],value*ss[:,0],'three-cell eigenpair residual',tol=2e-8)
        lower=len(edges)*np.linalg.eigvalsh(MN(.5,gg))[0];upper=np.linalg.eigvalsh(MN(1,gg))[0]
        le(lower,value,'Anderson lower bound');le(value,upper,'matching upper bound')
        brackets.append(dict(graph=name,g_over_Delta=gg,lower=float(lower),computed=value,upper=float(upper)))
# Negative controls make the scope failures concrete.
for false in [int(np.count_nonzero(np.abs(np.linalg.eigvalsh(H0N))<1e-10))==1,
              np.allclose(np.linalg.eigvalsh(MN(0,1))[1]-np.linalg.eigvalsh(MN(0,1))[0],np.linalg.eigvalsh(VN)[1]-np.linalg.eigvalsh(VN)[0]),
              s.Rational(1,8)==1]:
    if bool(false):raise RuntimeError('negative control accepted')
    counts['negative_controls']+=1
print(json.dumps(dict(status='passed',**counts,uniform_ratios=37,disordered_positive_link_cases=12,
    graph_cases=brackets,decoder_outcomes=len(branches),maximum_scaled_energy_residual=worst_energy,
    maximum_vector_infidelity_residual=worst_vector,minimum_disordered_coefficient_eigenvalue=disordered_min_rank_value,
    selected_uniform_cases=[reference[k] for k in [0,6,12,18,24,30,36]],
    python=platform.python_version(),sympy=s.__version__,numpy=np.__version__,scipy=scipy.__version__,
    scope='finite exact identities and bounded numerical tests; universal statements rely on the written proofs; no Lean, QPU, thermodynamic gap or historical-priority claim'),indent=2))
