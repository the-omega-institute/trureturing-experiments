"""Exact angular momentum and coherent vertex sampling for a unit-flux octahedron."""
import sympy as sp

half = sp.Rational(1,2)
eye = sp.eye
kron = sp.kronecker_product
root = sp.sqrt(2)
I = sp.I
pauli = (sp.Matrix([[0,1],[1,0]]), sp.Matrix([[0,-I],[I,0]]), sp.diag(1,-1))
sites = (sp.Matrix([1,0]),sp.Matrix([0,1]),sp.Matrix([1,1])/root,
         sp.Matrix([1,-1])/root,sp.Matrix([1,I])/root,sp.Matrix([1,-I])/root)


def spin_one():
    raise_op = sp.Matrix([[0,root,0],[0,0,root],[0,0,0]])
    return ((raise_op+raise_op.T)/2,(raise_op-raise_op.T)/(2*I),sp.diag(1,0,-1))


def simple(a):
    return a.applyfunc(sp.simplify)


rows = []
for z in sites:
    pair = sp.Matrix([z[0]**2,root*z[0]*z[1],z[1]**2])
    rows.append(kron(z.conjugate().T,pair.T))
S = sp.Matrix.vstack(*rows)
G = simple(S*S.H)
assert simple((G-2*eye(6))*(G-half*eye(6))) == sp.zeros(6)
U = simple(root*(7*eye(6)-2*G)*S/6)
assert simple(U*U.H) == eye(6)
l_ops = [kron(s/2,eye(3)) for s in pauli]
k_ops = [kron(eye(2),-s.T) for s in spin_one()]
g = (eye(12)+2*sum((kron(l,s) for l,s in zip(l_ops,pauli)),sp.zeros(12)))/2
h = -(eye(12)+2*sum((kron(k,s) for k,s in zip(k_ops,pauli)),sp.zeros(12)))/3
assert g*g == h*h == eye(12)
# Test the spectrum before sampling, keeping the entries in Q(sqrt(2),i).
base = simple(eye(12)-g*h)
assert simple(g*base+base*g-base*g*base) == sp.zeros(12)
assert (g+h).trace()/2 == 1
squared = simple(base.H*base)
t = sp.Symbol("t")
assert sp.factor(squared.charpoly(t).as_expr()) == sp.factor(t*(t-sp.Rational(4,3))**6*(t-4)**5)
# Every change of relative phase on the two irreducible summands has three
# Fourier blocks. Only the north to south entry is needed for locality.
P = simple((S.H*S-half*eye(6))*sp.Rational(2,3))
Q = eye(6)-P
assert simple(P*P-P)==sp.zeros(6)
np,nq,spole,sq = [kron(simple(U[i,:]*p),eye(2)) for i,p in ((0,P),(0,Q),(1,P),(1,Q))]
fixed = simple(np*base*spole.H + nq*base*sq.H)
forward,backward = simple(nq*base*spole.H),simple(np*base*sq.H)
ladder=sp.Matrix([[0,1],[0,0]])
assert fixed==ladder/3 and forward==-2*ladder/9 and backward==2*ladder/9
print("one zero mode, squared singular values 0, 4/3, 4 with multiplicities 1, 6, 5")
print("antipodal block = (1/3 - (4i/9) sin(theta)) times the spin raising matrix")
