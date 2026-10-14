#!/usr/bin/env python3
"""Exact finite HS4 exchange-readout and classical FIB closure checks.
No Lean or hardware claim. Run with Python and SymPy; no network needed.
"""
from __future__ import annotations
import itertools
import json
import sympy as s

checks=[]
def eq(label, a, b):
    d=a-b
    ok=all(s.simplify(v)==0 for v in d) if isinstance(d,s.MatrixBase) else s.simplify(d)==0
    if not ok:
        raise RuntimeError(label+': failed')
    checks.append(label)

def kron(*ms):
    out=s.Matrix([[1]])
    for m in ms: out=s.kronecker_product(out,m)
    return out

I=s.eye(2); X=s.Matrix([[0,1],[1,0]]); Y=s.Matrix([[0,-s.I],[s.I,0]]); Z=s.diag(1,-1)
paulis=[I,X,Y,Z]

def swap(i,j):
    m=s.zeros(16)
    for col,b in enumerate(itertools.product(range(2),repeat=4)):
        v=list(b);v[i],v[j]=v[j],v[i]
        row=sum(v[k]*2**(3-k) for k in range(4))
        m[row,col]=1
    return m

Sw={(i,j):swap(i,j) for i in range(4) for j in range(i+1,4)}
e0=s.zeros(16,1);e1=s.zeros(16,1)
e0[3]=e0[12]=1/s.sqrt(3)
for j in [5,6,9,10]: e0[j]=-1/(2*s.sqrt(3))
for j in [5,10]: e1[j]=s.Rational(1,2)
for j in [6,9]: e1[j]=-s.Rational(1,2)
E=e0.row_join(e1);P=E*E.H
hp=(e0+s.I*e1)/s.sqrt(2);hm=s.conjugate(hp)
eq('isometry',E.H*E,I)
eq('endpoint orthogonality',(hp.H*hm)[0],0)
for edge,m in Sw.items():eq('swap preserves singlets '+str(edge),m*E,E*(E.H*m*E))
eq('AB compression',E.H*Sw[0,1]*E,Z)
eq('AC compression',E.H*Sw[0,2]*E,-s.sqrt(3)*X/2-Z/2)
eq('AD compression',E.H*Sw[0,3]*E,s.sqrt(3)*X/2-Z/2)
cols=[]
for inds in itertools.product(range(4),repeat=4):
    if sum(i!=0 for i in inds)>2:continue
    M=kron(*(paulis[i] for i in inds));C=s.simplify(E.H*M*E)
    eq('two-local blind '+str(inds),s.trace(C*Y),0)
    cols.append(C.reshape(4,1))
eq('two-local restricted rank',s.Matrix.hstack(*cols).rank(),3)
A=(Sw[0,3]-Sw[0,2])/s.sqrt(3)
eq('exchange X compression',E.H*A*E,X)
eq('exchange cubic',A**3,A)
eq('exchange eigenvalues',s.Matrix(A).charpoly().as_expr(),s.Symbol('lambda')**8*(s.Symbol('lambda')-1)**4*(s.Symbol('lambda')+1)**4)
U=s.eye(16)+(1/s.sqrt(2)-1)*A**2-s.I*A/s.sqrt(2)
eq('one-pulse unitary',s.simplify(U.H*U),s.eye(16))
M=Sw[0,1]
eq('one pulse reads Y',s.simplify(E.H*U.H*M*U*E),Y)
eq('plus endpoint deterministic',s.simplify((hp.H*U.H*M*U*hp)[0]),1)
eq('minus endpoint deterministic',s.simplify((hm.H*U.H*M*U*hm)[0]),-1)
t=s.symbols('t',real=True)
u=s.cos(t)*I-s.I*s.sin(t)*X
eq('all-angle readout',s.trigsimp(u.H*Z*u),s.cos(2*t)*Z+s.sin(2*t)*Y)
eq('one commutator closure',s.I*(X*Z-Z*X),2*Y)
eq('dynamic rank',s.Matrix.hstack(I.reshape(4,1),X.reshape(4,1),Z.reshape(4,1),Y.reshape(4,1)).rank(),4)
K=2*s.I*(Sw[0,1]*Sw[1,2]-Sw[1,2]*Sw[0,1])
eq('chirality compression',s.simplify(E.H*K*E),-2*s.sqrt(3)*Y)
eq('chirality cubic',K**3,12*K)
for sign in (-1,1):
    Us=s.eye(16)+(1/s.sqrt(2)-1)*A**2-sign*s.I*A/s.sqrt(2)
    eq('pulse orientation '+str(sign),s.simplify(E.H*Us.H*M*Us*E),sign*Y)
# Single-exchange alternative, with contrast bounded by sqrt(3)/2.
w=s.cos(t)*I-s.I*s.sin(t)*Z
C=E.H*Sw[0,2]*E
single=s.trigsimp(w.H*C*w)
eq('one elementary exchange contrast',s.simplify(s.trace(Y*single)/2),s.sqrt(3)*s.sin(2*t)/2)
# Measurement-only histories: real projectors and swap gates remain blind.
for n in range(1,5):
    for labels in itertools.product((0,1),repeat=n):
        L=I
        for k,o in enumerate(labels):
            S=E.H*Sw[(0,1) if k%2==0 else (0,2)]*E
            L=(I+(-1)**o*S)/2*L
        eq('record-only blind '+str(labels),s.trace(Y*L.H*L),0)
# Classical observation quotient: order is 0,1,2,3,13; coords x=site1,y=site3,z=site2.
ones=s.ones(5,1);x=s.Matrix([0,1,0,0,1]);y=s.Matrix([0,0,0,1,1]);z=s.Matrix([0,0,1,0,0]);d=s.Matrix([1,-1,0,-1,1]);xy=x.multiply_elementwise(y)
B=s.Matrix.hstack(ones,x,y,z)
eq('FIB static rank',B.rank(),4)
eq('FIB blind direction',B.T*d,s.zeros(4,1))
eq('FIB product reads hidden',(d.T*xy)[0],1)
eq('FIB completed rank',B.row_join(xy).rank(),5)
a,b,c,e,aa,bb,cc,ee=s.symbols('a b c e aa bb cc ee',real=True)
f=a*ones+b*x+c*y+e*z;g=aa*ones+bb*x+cc*y+ee*z
eq('FIB symmetric multiplication defect',(d.T*f.multiply_elementwise(g))[0],b*cc+c*bb)
F=a*I+b*X+c*Z;G=aa*I+bb*X+cc*Z
eq('HS Lie defect coefficient',s.trace(Y*s.I*(F*G-G*F))/2,2*(b*cc-c*bb))
eq('HS Jordan closure',s.trace(Y*(F*G+G*F))/2,0)
r=s.symbols('r',real=True);rho=(I+r*Y)/2
q=(1+r*s.sin(2*t))/2
fisher=s.diff(q,r)**2/q+s.diff(1-q,r)**2/(1-q)
eq('Fisher response',fisher,s.sin(2*t)**2/(1-r*r*s.sin(2*t)**2))
eq('optimal outcome probability',s.trace(rho*(I+Y)/2),(1+r)/2)
# Tetrahedral separable pair realization and quantitative gluing obstruction.
def ptrace(v, keep):
    out=s.zeros(2**len(keep)); other=[k for k in range(4) if k not in keep]
    for bits in itertools.product(range(2),repeat=len(other)):
        vv=[]
        for kept in itertools.product(range(2),repeat=len(keep)):
            b=[0]*4
            for k,a in zip(keep,kept):b[k]=a
            for k,a in zip(other,bits):b[k]=a
            vv.append(v[sum(b[k]*2**(3-k) for k in range(4))])
        w=s.Matrix(vv);out+=w*w.H
    return s.simplify(out)
ss=s.Matrix([0,1,-1,0])/s.sqrt(2);R=s.eye(4)/6+ss*ss.H/3
for edge in Sw:
    eq('HS pair R '+str(edge),ptrace(hp,edge),R)
verts=[s.Matrix(v)/s.sqrt(3) for v in [(1,1,1),(1,-1,-1),(-1,1,-1),(-1,-1,1)]]
def bloch(v):return (I+sum((v[j]*paulis[j+1] for j in range(3)),s.zeros(2)))/2
sep_pair=s.zeros(4)
for n in verts:
    eq('tetra norm '+str(tuple(n)),(n.T*n)[0],1)
    sep_pair+=kron(bloch(n),bloch(-n))/4
eq('tetra separable pair',s.simplify(sep_pair),R)
Dot=kron(X,X)+kron(Y,Y)+kron(Z,Z)
Rp=(s.eye(4)-Dot/9)/4
eq('closest separable pair spectrum',s.factor((R-Rp).charpoly().as_expr()),(s.Symbol('lambda')-s.Rational(1,6))*(s.Symbol('lambda')+s.Rational(1,18))**3)
patterns=[sg for sg in itertools.product((-1,1),repeat=4) if sum(sg)==0]
for i,j in Sw:
    eq('balanced pair moment '+str((i,j)),sum(sg[i]*sg[j] for sg in patterns)/s.Integer(6),-s.Rational(1,3))
eq('tetra second moment',sum((n*n.T for n in verts),s.zeros(3))/4,s.eye(3)/3)
eq('separable trace norm lower bound is attained',6*(s.Rational(1,6)+3*s.Rational(1,18)),2)
# The exact oriented triple tensor fixes all projective three-party statistics.
for a,b,c in itertools.product(range(3),repeat=3):
    M3=kron(paulis[a+1],paulis[b+1],paulis[c+1],I)
    eq('oriented triple '+str((a,b,c)),s.simplify((hp.H*M3*hp)[0]),-s.LeviCivita(a,b,c)/s.sqrt(3))
# Independent fair edge signs, and one private local sign, reproduce orthogonal-axis outcomes.
v=s.symbols('v',real=True)
for xx,yy,zz in itertools.product((-1,1),repeat=3):
    mass=s.Integer(0)
    for al,be,ga,eta in itertools.product((-1,1),repeat=4):
        if (be*ga,ga*al,eta*al*be)==(xx,yy,zz):mass+=(1+eta*v)/16
    eq('triangle local parity '+str((xx,yy,zz)),mass,(1+v*xx*yy*zz)/8)
# Discarding the signed pulse label destroys the exact hidden signal.
for o in (-1,1):
    eq('discarded pulse sign '+str(o),sum((1+sg*o*r*s.sin(2*t))/4 for sg in (-1,1)),s.Rational(1,2))
# Negative controls use explicit false formulas, separately from count above.
negatives=[s.simplify(s.trace(rho*Z)-r),s.simplify(s.trace(E.H*U.H*M*U*E*Y)/2+1),s.Integer(3-4)]
if not all(v!=0 for v in negatives):raise RuntimeError('negative control failed')
print(json.dumps({'status':'passed','exact_checks':len(checks),'negative_controls':len(negatives),'static_observation_rank':3,'one_pulse_observation_rank':4,'fib_static_rank':4,'fib_product_rank':5,'pulse_angle':'pi/4','pulse_probability_plus':'(1+r)/2','single_swap_success':'(2+sqrt(3))/4','sympy':s.__version__,'scope':'exact finite and symbolic identities; no Lean compilation; no all-history proof by enumeration'},indent=2))
