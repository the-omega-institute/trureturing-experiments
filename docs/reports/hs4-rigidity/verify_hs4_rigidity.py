#!/usr/bin/env python3
"""Exact HS4 rigidity/parent-Hamiltonian checks and bounded numerical tests.

This script does not recompile Lean or certify the global orbit-classification proof.
Exact checks use SymPy; numerical tests use a fixed NumPy seed. Run without arguments.
"""
import itertools
import json
import math
import sys
import numpy as np
import sympy as sp

checks = []

def require(condition, label):
    if not bool(condition):
        raise AssertionError(label)
    checks.append(label)

def zero(matrix):
    return all(sp.simplify(x) == 0 for x in matrix)

def kron(mats):
    value = sp.Matrix([[1]])
    for matrix in mats:
        value = sp.kronecker_product(value, matrix)
    return value

def permutation(order):
    value = sp.zeros(16)
    for index in range(16):
        bits = [(index >> (3-j)) & 1 for j in range(4)]
        target = sum(bits[order[j]] << (3-j) for j in range(4))
        value[target, index] = 1
    return value

def flatten(state, keep):
    drop = [j for j in range(4) if j not in keep]
    value = sp.zeros(2**len(keep), 2**len(drop))
    for index in range(16):
        bits = [(index >> (3-j)) & 1 for j in range(4)]
        r = sum(bits[k] << (len(keep)-1-j) for j, k in enumerate(keep))
        c = sum(bits[k] << (len(drop)-1-j) for j, k in enumerate(drop))
        value[r, c] = state[index]
    return value

R, I = sp.Rational, sp.I
identity = sp.eye(2)
pauli = [sp.Matrix([[0, 1], [1, 0]]), sp.Matrix([[0, -I], [I, 0]]), sp.diag(1, -1)]
omega = (-1 + I*sp.sqrt(3))/2
plus = sp.zeros(16, 1)
for index, coefficient in {3:1, 12:1, 5:omega, 10:omega, 6:sp.conjugate(omega), 9:sp.conjugate(omega)}.items():
    plus[index] = coefficient/sp.sqrt(6)
minus = plus.conjugate()
require(zero(plus.conjugate().T*plus-sp.ones(1)), 'HS normalized')
require(zero(plus.conjugate().T*minus), 'conjugates orthogonal')
pair = sp.Matrix([[R(1,6),0,0,0], [0,R(1,3),-R(1,6),0], [0,-R(1,6),R(1,3),0], [0,0,0,R(1,6)]])
t = sp.Symbol('t', real=True)
for keep in itertools.combinations(range(4), 2):
    matrices = [flatten(state, keep) for state in (plus, minus)]
    for sign, matrix in zip(('+','-'), matrices):
        require(zero(matrix*matrix.conjugate().T-pair), 'pair '+str(keep)+' '+sign)
    require(zero(t*matrices[0]*matrices[0].conjugate().T+(1-t)*matrices[1]*matrices[1].conjugate().T-pair), 'mixed fiber '+str(keep))

cut_matrices = [flatten(plus, keep) for keep in ((0,1),(0,2),(0,3))]
determinants = [sp.simplify(matrix.det()) for matrix in cut_matrices]
expected = [I*sp.sqrt(3)/36, R(1,24)+I*sp.sqrt(3)/72, R(1,24)-I*sp.sqrt(3)/72]
require(all(sp.simplify(a-b)==0 for a,b in zip(determinants, expected)), 'literal cut determinants')
chi = sp.simplify(sp.im(determinants[0]*sp.conjugate(determinants[1])))
require(chi == sp.sqrt(3)/864, 'LU invariant chirality value')
require(zero(permutation((0,1,3,2))*plus-minus), 'C D swap exchanges conjugates')

e1, e0 = sp.zeros(16,1), sp.zeros(16,1)
for k,v in {5:1,6:-1,9:-1,10:1}.items(): e1[k]=R(v,2)
for k,v in {3:1,12:1,5:-R(1,2),6:-R(1,2),9:-R(1,2),10:-R(1,2)}.items(): e0[k]=v/sp.sqrt(3)
E = sp.Matrix.hstack(e0,e1)
require(zero(E.T*E-sp.eye(2)), 'real singlet basis orthonormal')
require(zero(plus-(e0+I*e1)/sp.sqrt(2)), 'HS singlet coordinates')
J = [sum((kron([pauli[a] if j==q else identity for j in range(4)]) for q in range(4)), sp.zeros(16))/2 for a in range(3)]
L = sum((j*j for j in J), sp.zeros(16))
require(zero(L*E), 'singlet basis in kernel')
require(L.rank()==14, 'total spin kernel dimension two')
require(zero(L*(L-2*sp.eye(16))*(L-6*sp.eye(16))), 'total spin polynomial')
AB, BC, AC = permutation((1,0,2,3)), permutation((0,2,1,3)), permutation((2,1,0,3))
require(zero(E.T*AB*E-sp.diag(1,-1)), 'AB swap compression')
require(zero(E.T*AC*E-sp.Matrix([[-R(1,2),-sp.sqrt(3)/2],[-sp.sqrt(3)/2,R(1,2)]])), 'AC swap compression')
K = sp.zeros(16)
for p in itertools.permutations(range(3)):
    sign = (-1)**sum(p[a]>p[b] for a in range(3) for b in range(a+1,3))
    K += sign*kron([pauli[p[0]],pauli[p[1]],pauli[p[2]],identity])
require(zero(K-K.conjugate().T), 'scalar chirality Hermitian')
require(zero(K-2*I*(AB*BC-BC*AB)), 'chirality permutation identity')
require(zero(K**3-12*K), 'chirality polynomial')
require(zero(K*plus+2*sp.sqrt(3)*plus), 'HS plus chirality eigenvalue')
require(zero(K*minus-2*sp.sqrt(3)*minus), 'HS minus chirality eigenvalue')
H = (L+sp.eye(16)+K/(2*sp.sqrt(3))).applyfunc(sp.simplify)
require(zero(H*plus), 'parent annihilates HS plus')
require(H.eigenvals()=={0:1,2:4,3:3,4:3,7:5}, 'parent exact spectrum and gap')

bell = [sp.Matrix([1,0,0,1])/sp.sqrt(2), sp.Matrix([1,0,0,-1])/sp.sqrt(2), sp.Matrix([0,1,1,0])/sp.sqrt(2), sp.Matrix([0,1,-1,0])/sp.sqrt(2)]
U = sp.Matrix.hstack(*(sp.kronecker_product(b,b) for b in bell))
perms = [(2,3,0,1),(1,0,3,2)]
for operator, p in zip([kron([pauli[0],identity,pauli[0],identity]),kron([pauli[2],identity,pauli[2],identity])],perms):
    target = sp.zeros(4)
    for j in range(4): target[p[j],j]=1
    require(zero(operator*U-U*target), 'Bell-coordinate permutation '+str(p))
for operator, signs in [(kron([pauli[0],pauli[0],identity,identity]),(1,-1,1,-1)),(kron([pauli[2],pauli[2],identity,identity]),(1,1,-1,-1))]:
    require(zero(operator*U-U*sp.diag(*signs)), 'Bell-coordinate signs '+str(signs))

x = sp.Symbol('x', positive=True)
a,b = R(1,6),R(1,2)
p = (4-3*sp.log(3))/8+((-11+2*sp.log(2)+12*sp.log(3))/2)*x+((36-39*sp.log(3))/2)*x*x+18*(sp.log(3)-1)*x**3
g = p+x*sp.log(x)
for point in (a,b):
    require(sp.simplify(g.subs(x,point))==0, 'entropy contact '+str(point))
    require(sp.simplify(sp.diff(g,x).subs(x,point))==0, 'entropy derivative contact '+str(point))
require(sp.simplify(sp.diff(g,x,4)-2/x**3)==0, 'gap fourth derivative')
A=(11+2*sp.log(2)-9*sp.log(3))/2
B=3*(5*sp.log(3)-6)/2
C=54*(sp.log(3)-1)
require(sp.simplify(3*A+B+C*R(5,36)-3*sp.log(2)-R(3,2)*sp.log(3))==0, 'gap constant identity')
require(9*12**3<25**3 and 37**9>32*25**9, 'Renyi 5/3 rational comparison')
require(3*15**2<28**2, 'Renyi 3/2 square comparison')

# Negative controls are genuine rejected exact identities.
negative_controls = {
    'wrong chirality sign': zero(K*plus-2*sp.sqrt(3)*plus),
    'collapsed LU invariant branches': sp.simplify(chi+chi)==0,
    'unconjugated Gram': zero(cut_matrices[0]*cut_matrices[0].T-pair),
}
require(not any(negative_controls.values()), 'all negative controls rejected')

rng=np.random.default_rng(20261014)
bound=1+0.5*math.log2(3)
Af,Bf,Cf=[float(v) for v in (A,B,C)]
pfun=sp.lambdify(x,p,'numpy')
states=[]
for _ in range(512):
    z=rng.normal(size=16)+1j*rng.normal(size=16);states.append(z/np.linalg.norm(z))
base=np.array(plus.evalf(),dtype=complex).ravel()
for exponent in (1,2,3,4,5,6):
    for _ in range(8):
        z=base+10**(-exponent)*(rng.normal(size=16)+1j*rng.normal(size=16));states.append(z/np.linalg.norm(z))
spinflip=np.array(kron([pauli[1]]*4),dtype=complex)
max_gap_residual=0.0
for index,z in enumerate(states):
    spectra=[]
    for keep in ((0,1),(0,2),(0,3)):
        drop=tuple(j for j in range(4) if j not in keep)
        f=z.reshape((2,)*4).transpose(keep+drop).reshape(4,4)
        spectra.append(np.linalg.eigvalsh(f@f.conj().T))
    spectra=np.maximum(np.array(spectra),0)
    p2=float(np.sum(spectra*spectra))
    e3=sum(float(np.prod(v[list(sub)])) for v in spectra for sub in itertools.combinations(range(4),3))
    logs=np.zeros_like(spectra);mask=spectra>0;logs[mask]=spectra[mask]*np.log(spectra[mask])
    deficit=bound+float(np.sum(logs))/(3*math.log(2))
    rhs=float(np.sum(pfun(spectra)+logs))-Bf*(p2-1)+Cf*(5/36-e3)
    residual=abs(3*math.log(2)*deficit-rhs);max_gap_residual=max(max_gap_residual,residual)
    require(residual<1e-11, 'numerical gap identity '+str(index))
    require(p2>=1-1e-11 and e3<=5/36+1e-11 and deficit>=-1e-11, 'numerical inherited bounds '+str(index))
    spectral_residual=float(np.sum((spectra-1/6)**2*(spectra-1/2)**2))
    require(spectral_residual<=36*math.log(2)*deficit+1e-10, 'numerical spectral stability '+str(index))
    bloch_sq=0.0
    for q in range(4):
        f=z.reshape((2,)*4).transpose((q,)+tuple(j for j in range(4) if j!=q)).reshape(2,8)
        rho=f@f.conj().T
        bloch_sq+=2*float(np.trace(rho@rho).real)-1
    tangle=abs(z.T@spinflip@z)**2
    require(abs(2*(p2-1)-bloch_sq-tangle)<1e-11, 'numerical purity decomposition '+str(index))
print(json.dumps({'status':'passed','exact_checks':len(checks)-4*len(states),'numerical_states':len(states),'numerical_checks':4*len(states),'negative_controls':len(negative_controls),'maximum_gap_identity_residual':max_gap_residual,'parent_spectrum':{'0':1,'2':4,'3':3,'4':3,'7':5},'chirality_invariant':'sqrt(3)/864','sympy':sp.__version__,'numpy':np.__version__,'python':sys.version.split()[0],'scope':'exact finite identities plus bounded numerical tests; no Lean replay; no universal-proof certification'},indent=2))
