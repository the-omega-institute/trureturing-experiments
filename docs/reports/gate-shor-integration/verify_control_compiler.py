#!/usr/bin/env python3
"""Native no11 controls: exact identities, modular rank certificates, pulse compiler.
No hardware or Lean execution. All finite checks have explicit scope.
"""
from __future__ import annotations
import collections
import hashlib
import itertools
import json
import math
import platform
import numpy as np
import sympy as sp

COUNTS = {'exact_assertions': 0, 'modular_rank_assertions': 0,
          'numerical_assertions': 0, 'negative_controls': 0}
def require(ok, label, category='exact_assertions'):
    if not bool(ok):
        raise RuntimeError(label)
    COUNTS[category] += 1

def exact_equal(a, b, label):
    delta = a-b
    vals = list(delta) if isinstance(delta, sp.MatrixBase) else [delta]
    require(all(sp.simplify(x) == 0 for x in vals), label)

def close(a, b, label, tol=2e-10):
    require(np.max(np.abs(np.asarray(a)-np.asarray(b))) <= tol,
            label, 'numerical_assertions')

def legal_words(n):
    return [w for w in itertools.product((0,1), repeat=n)
            if all(w[k]*w[k+1] == 0 for k in range(n-1))]

def toggle_data(n):
    words = legal_words(n); idx = {w:j for j,w in enumerate(words)}
    edge_sets = []
    for k in range(n):
        edges = []
        for a,w in enumerate(words):
            if w[k] == 0 and (k == 0 or w[k-1] == 0) and (k == n-1 or w[k+1] == 0):
                v=list(w);v[k]=1;edges.append((a,idx[tuple(v)]))
        edge_sets.append(edges)
    return words, edge_sets

def generators(n):
    words,edges=toggle_data(n);d=len(words);xs=[];ks=[]
    for es in edges:
        x=np.zeros((d,d),dtype=np.int64);k=x.copy()
        for a,b in es:x[a,b]=x[b,a]=1;k[a,b]=-1;k[b,a]=1
        xs.append(x);ks.append(k)  # Y = i K; complexified Lie generators can use X,K
    return words,edges,xs,ks

# Incremental finite-field elimination. Every stored raw matrix is an actual
# nested commutator of the original integer generators, never a numerical fit.
def lie_rank(n, phase_site=None, prime=1009):
    words,edges,xs,ks=generators(n);d=len(words);gs=xs+([] if phase_site is None else [ks[phase_site]])
    basis={};raw=[];ancestry=[]
    def add(mat,parent):
        vec=mat.reshape(-1).copy()%prime
        for p,row in sorted(basis.items()):
            if vec[p]:vec=(vec-int(vec[p])*row)%prime
        nz=np.flatnonzero(vec)
        if len(nz)==0:return False
        p=int(nz[0]);basis[p]=(vec*pow(int(vec[p]),-1,prime))%prime
        raw.append(mat.copy()%prime);ancestry.append(parent);return True
    for j,g in enumerate(gs):add(g,[-1,j])
    q=0
    while q<len(raw) and len(raw)<d*d-1:
        for j,g in enumerate(gs):
            add((g@raw[q]-raw[q]@g)%prime,[q,j])
            if len(raw)==d*d-1:break
        q+=1
    # Independent ancestry replay: raw[i] is a Lie word over integer matrices.
    for i,(parent,j) in enumerate(ancestry):
        expect=gs[j]%prime if parent<0 else (gs[j]@raw[parent]-raw[parent]@gs[j])%prime
        require(np.array_equal(expect,raw[i]),'commutator ancestry','modular_rank_assertions')
    require(all(int(np.trace(a))%prime==0 for a in raw),'trace-zero upper bound','modular_rank_assertions')
    if phase_site is not None:
        require(len(raw)==d*d-1,'full complexified Lie rank','modular_rank_assertions')
    else:
        require(len(raw)<=d*(d-1)//2,'orthogonal-form upper bound','modular_rank_assertions')
    digest=hashlib.sha256(json.dumps(ancestry,separators=(',',':')).encode()).hexdigest()
    return {'n':n,'dimension':d,'phase_site_one_based':None if phase_site is None else phase_site+1,
            'rank_mod_prime':len(raw),'prime':prime,'word_ancestry_sha256':digest}

I=sp.eye(2);X=sp.Matrix([[0,1],[1,0]]);Y=sp.Matrix([[0,-sp.I],[sp.I,0]]);Z=sp.diag(1,-1)
def symbolic_controls(n):
    words,edges,xs,ks=generators(n)
    return words,edges,[sp.Matrix(x) for x in xs],[sp.I*sp.Matrix(k) for k in ks]

def symbolic_pulse(h,t):
    return sp.eye(h.rows)+(sp.cos(t)-1)*(h*h)-sp.I*sp.sin(t)*h

def code_isometry(m):
    n=2*m-1;words=legal_words(n);idx={w:j for j,w in enumerate(words)}
    e=sp.zeros(len(words),2**m)
    for a,bits in enumerate(itertools.product((0,1),repeat=m)):
        w=tuple(x for j,b in enumerate(bits) for x in ([b,0] if j<m-1 else [b]))
        e[idx[w],a]=1
    return e

def qubit_gate_matrix(m,q,op):
    return sp.kronecker_product(*[op if j==q else I for j in range(m)])

def cnot_matrix(m,c,t):
    ans=sp.zeros(2**m)
    for x in range(2**m):
        y=x^((1<<(m-1-t)) if (x>>(m-1-c))&1 else 0);ans[y,x]=1
    return ans

# Symbolic all-input identities on m=2 and m=3 encoded registers.
for m in [2,3]:
    n=2*m-1;words,edges,hx,hy=symbolic_controls(n);d=len(words);E=code_isometry(m)
    exact_equal(E.T*E,sp.eye(2**m),'encoding isometry')
    angle=sp.symbols('t',real=True)
    for q in range(m):
        for h,p in [(hx[2*q],X),(hy[2*q],Y)]:
            native=symbolic_pulse(h,angle)
            target=sp.cos(angle)*I-sp.I*sp.sin(angle)*p
            exact_equal(native*E,E*qubit_gate_matrix(m,q,target),'symbolic single-qubit intertwiner')
    for q in range(m-1):
        mid=2*q+1
        D=symbolic_pulse(hx[mid],sp.pi)
        wanted=sp.eye(2**m)
        for a in range(2**m):
            if ((a>>(m-1-q))&1)==0 and ((a>>(m-2-q))&1)==0:wanted[a,a]=-1
        exact_equal(D*E,E*wanted,'conditional minus phase and clean separator')
        # S is applied first, then D, then S^dagger.
        S=symbolic_pulse(hx[2*q+2],sp.pi/2)*symbolic_pulse(hx[2*q],sp.pi/2)
        CZ=S.H*D*S
        cz=sp.eye(2**m)
        for a in range(2**m):
            if ((a>>(m-1-q))&1) and ((a>>(m-2-q))&1):cz[a,a]=-1
        exact_equal(CZ*E,E*cz,'five-pulse exact CZ')
        for c,t in [(q,q+1),(q+1,q)]:
            h=symbolic_pulse(hx[2*t],sp.pi/2)*symbolic_pulse(hy[2*t],sp.pi/4)
            cx=h*CZ*h.H
            exact_equal(cx*E,E*cnot_matrix(m,c,t),'nine-pulse exact CNOT')

# The phase-zero menu has a conjugated-real orthogonal representation.
for n in range(1,6):
    words,edges,hx,hy=symbolic_controls(n);d=len(words)
    W=sp.diag(*[sp.I**sum(w) for w in words]);parity=sp.diag(*[(-1)**sum(w) for w in words])
    for h in hx:
        a=sp.simplify(W.H*(-sp.I*h)*W)
        exact_equal(a,sp.conjugate(a),'real gauged generator')
        exact_equal(a.T,-a,'antisymmetric gauged generator')
        exact_equal(parity*h*parity,-h,'parity flips phase-zero Hamiltonian')
    u=words.index((0,)*n);v=words.index((1,)+(0,)*(n-1))
    plus=sp.zeros(d,1);minus=sp.zeros(d,1)
    plus[u]=plus[v]=1/sp.sqrt(2);minus[u]=1/sp.sqrt(2);minus[v]=-1/sp.sqrt(2)
    exact_equal(parity*sp.conjugate(plus),minus,'antiunitary blind pair')
    exact_equal((plus.H*minus)[0],0,'orthogonal blind hypotheses')
    test=symbolic_pulse(hy[0],sp.pi/4)
    exact_equal((test*plus)[v],1,'one quadrature gives deterministic first outcome')
    exact_equal((test*minus)[u],1,'one quadrature gives deterministic second outcome')

class NativeCompiler:
    """Known H/P/CNOT circuit -> specified local pulses on 2m-1 no11 sites.
    Phase records a known constant relating the native and named gates.
    """
    def __init__(self,m):
        self.m=m;self.n=2*m-1;self.words,self.edges=toggle_data(self.n)
        self.E=np.array(code_isometry(m),dtype=complex);self.state=self.E.copy()
        self.phase=1+0j;self.pulses=[];self.native_cnot_count=0
    def pulse(self,site,beta,t):
        old=self.state.copy();c=np.cos(t);s=np.sin(t)
        for a,b in self.edges[site]:
            self.state[a]=c*old[a]-1j*s*np.exp(-1j*beta)*old[b]
            self.state[b]=c*old[b]-1j*s*np.exp(1j*beta)*old[a]
        self.pulses.append([site+1,float(beta),float(t)])
    def rx(self,q,t):self.pulse(2*q,0,t)
    def ry(self,q,t):self.pulse(2*q,np.pi/2,t)
    def rz(self,q,t):
        self.rx(q,-np.pi/4);self.ry(q,t);self.rx(q,np.pi/4)
    def h(self,q):
        self.ry(q,np.pi/4);self.rx(q,np.pi/2);self.phase*=-1j
    def p(self,q,t):
        self.rz(q,t/2);self.phase*=np.exp(-1j*t/2)
    def local_cnot(self,c,t):
        require(abs(c-t)==1,'adjacent CNOT routing')
        # h^dagger on target, S, D, S^dagger, h on target.
        self.rx(t,-np.pi/2);self.ry(t,-np.pi/4)
        q=min(c,t);self.rx(q,np.pi/2);self.rx(q+1,np.pi/2)
        self.pulse(2*q+1,0,np.pi)
        self.rx(q+1,-np.pi/2);self.rx(q,-np.pi/2)
        self.ry(t,np.pi/4);self.rx(t,np.pi/2)
        self.native_cnot_count+=1
    def swap_adj(self,a,b):
        self.local_cnot(a,b);self.local_cnot(b,a);self.local_cnot(a,b)
    def cx(self,c,t):
        if c==t:raise ValueError('distinct CNOT bits required')
        step=1 if t>c else -1;swaps=[];pos=c
        while abs(pos-t)>1:
            nxt=pos+step;self.swap_adj(pos,nxt);swaps.append((pos,nxt));pos=nxt
        self.local_cnot(pos,t)
        for a,b in reversed(swaps):self.swap_adj(a,b)
    def execute(self,ops):
        for op in ops:
            if op[0]=='H':self.h(op[1])
            elif op[0]=='P':self.p(op[1],op[2])
            elif op[0]=='CX':self.cx(op[1],op[2])
            else:raise ValueError(op)
        return self.state

def source_unitary(m,ops):
    u=np.eye(2**m,dtype=complex)
    for op in ops:
        if op[0]=='H':g=np.array(qubit_gate_matrix(m,op[1],(X+Z)/sp.sqrt(2)),dtype=complex)
        elif op[0]=='P':g=np.array(qubit_gate_matrix(m,op[1],sp.diag(1,sp.exp(sp.I*op[2]))).evalf(),dtype=complex)
        elif op[0]=='CX':g=np.array(cnot_matrix(m,op[1],op[2]),dtype=complex)
        else:raise ValueError(op)
        u=g@u
    return u

def ccx_ops(a,b,t):
    p=np.pi/4
    return [('H',t),('CX',b,t),('P',t,-p),('CX',a,t),('P',t,p),('CX',b,t),
            ('P',t,-p),('CX',a,t),('P',b,p),('P',t,p),('H',t),('CX',a,b),
            ('P',a,p),('P',b,-p),('CX',a,b)]
def cp_ops(a,b,theta):
    return [('P',a,theta/2),('P',b,theta/2),('CX',a,b),('P',b,-theta/2),('CX',a,b)]
def swap_ops(a,b):return [('CX',a,b),('CX',b,a),('CX',a,b)]
def inverse_ops(ops):
    return [('P',op[1],-op[2]) if op[0]=='P' else op for op in reversed(ops)]

# Independent Toffoli check against the exact permutation.
toff=source_unitary(3,ccx_ops(0,1,2));wanted=np.eye(8);wanted[:,[6,7]]=wanted[:,[7,6]]
close(toff,wanted,'standard Toffoli decomposition')
# Five logical qubits: 3-bit clock and the two-bit Fibonacci map mod 2.
ops=[('H',q) for q in range(3)]
for q in range(3):
    exponent=2**(2-q)
    for _ in range(exponent):
        ops+=ccx_ops(q,4,3)+ccx_ops(q,3,4)  # (x,y)->(x+y,x) mod 2
forward=[]
for q in range(3):
    forward.append(('H',q))
    for k in range(q+1,3):forward+=cp_ops(k,q,np.pi/(2**(k-q)))
forward+=swap_ops(0,2)
ops+=inverse_ops(forward)
logical=source_unitary(5,ops)
# Direct clock-and-permutation matrix; no gate decomposition reused here.
U=np.zeros((4,4),complex)
for x in range(2):
    for y in range(2):U[2*((x+y)%2)+x,2*x+y]=1
H=np.array([[1,1],[1,-1]],complex)/np.sqrt(2)
Hc=np.kron(np.kron(H,H),H)
control=np.zeros((32,32),complex)
for a in range(8):control[4*a:4*a+4,4*a:4*a+4]=np.linalg.matrix_power(U,a)
Fi=np.exp(-2j*np.pi*np.outer(np.arange(8),np.arange(8))/8)/np.sqrt(8)
direct=np.kron(Fi,np.eye(4))@control@np.kron(Hc,np.eye(4))
close(logical,direct,'complete source phase-estimation unitary')
compiler=NativeCompiler(5);compiled=compiler.execute(ops)
residual=float(np.linalg.norm(compiled-compiler.phase*compiler.E@direct,2))
close(compiled,compiler.phase*compiler.E@direct,'all-input native phase-estimation intertwiner')
close(compiled.conj().T@compiled,np.eye(32),'native circuit isometry')
leak=float(np.linalg.norm(compiled-compiler.E@(compiler.E.conj().T@compiled),2))
close(leak,0,'all separator bits reset at circuit end')
column=compiled[:,2];physical_prob=np.abs(column)**2;clock_prob=np.zeros(8)
for w,p in zip(compiler.words,physical_prob):
    codebits=w[::2];clock=4*codebits[0]+2*codebits[1]+codebits[2];clock_prob[clock]+=p
ideal_prob=np.sum(np.abs(direct[:,2].reshape(8,4))**2,axis=1)
close(clock_prob,ideal_prob,'compiled finite-clock distribution')

# Negative controls: two orthogonal states remain blind with one quadrature;
# a pi/2 separator pulse leaks out of the clean encoding; a fixed spectral
# readout cannot resolve coherence between different eigenvectors.
_,_,blind_x,_=symbolic_controls(3)
words3=legal_words(3);left=np.zeros(5,complex);right=left.copy()
i0=words3.index((0,0,0));i1=words3.index((1,0,0))
left[i0]=left[i1]=1/np.sqrt(2);right[i0]=1/np.sqrt(2);right[i1]=-1/np.sqrt(2)
readout=np.eye(5,dtype=complex)
for h,t in zip(blind_x,[sp.pi/4,sp.pi/7,sp.pi/5]):
    readout=np.array(symbolic_pulse(h,t).evalf(),dtype=complex)@readout
require(np.allclose(np.abs(readout@left)**2,np.abs(readout@right)**2,atol=2e-12),
        'phase-zero history fails to distinguish orthogonal pair','negative_controls')
_,_,xx,yy=symbolic_controls(3);Ee=code_isometry(2)
half=symbolic_pulse(xx[1],sp.pi/2)*Ee
require(half!=Ee*(Ee.T*half),'half pulse leaves separator record','negative_controls')
rplus=np.array([[.5,.5],[.5,.5]]);rminus=np.array([[.5,-.5],[-.5,.5]])
require(np.array_equal(np.diag(rplus),np.diag(rminus)) and not np.array_equal(rplus,rminus),
        'equal spectral weights do not identify the state','negative_controls')

lie_results=[]
for n in range(1,6):
    lie_results.append(lie_rank(n,None));lie_results.append(lie_rank(n,0))
for seed in [1,2]:lie_results.append(lie_rank(3,seed))

result={'status':'passed',**COUNTS,'rank_certificates':lie_results,
        'compiled_order_probe':{'ordinary_qubits':5,'native_sites':9,'legal_hilbert_dimension':len(compiler.words),
        'clock_size':8,'target_map':'(x,y)->(x+y,x) mod 2','input_target':[1,0],
        'source_H_P_CX_count':len(ops),'routed_adjacent_CNOTs':compiler.native_cnot_count,
        'native_pulses':len(compiler.pulses),'all_input_operator_residual':residual,'final_leakage_norm':leak,
        'clock_probabilities':clock_prob.tolist(),'uses_repeated_small_powers':True},
        'python':platform.python_version(),'numpy':np.__version__,'sympy':sp.__version__,
        'scope':'symbolic finite identities; exact finite-field Lie-word lower bounds; floating ideal-pulse simulation; no Lean, hardware, large-N factoring or speedup certificate'}
print(json.dumps(result,indent=2))
