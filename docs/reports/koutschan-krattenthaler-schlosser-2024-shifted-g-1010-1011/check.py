"""Exact checks for Koutschan-Krattenthaler-Schlosser, arXiv:2401.08481v1, Conjecture 23, identities (10.10)-(10.12):
the uniform column-difference identity (alpha,beta,gamma,delta) -> (alpha+1,beta+2,gamma,delta-4) under h, the shared-column
block structure of each link (common columns j >= m, zero upper-right block, leading-block determinants), the Proposition 22
pieces, and the identities themselves for n = 3..10."""
from fractions import Fraction as F
import sympy as sp
def binom(x,p):
    if p<0: return F(0)
    r=F(1)
    for s in range(p): r*=F(x-s)
    for s in range(1,p+1): r/=s
    return r
def entry(P,i,j):
    a,b,g,d=P
    return F(2)**(i+b)*binom(i+4*j+g,4*j+a)+binom(-i+4*j+d,4*j+a)
def A(P,N): return sp.Matrix(N,N,lambda i,j: sp.Rational(entry(P,i,j)))
def G(P,N): return F(1) if N<=0 else F(str(A(P,N).det()))
t=sp.symbols('t')
h=1/((1-t)*(1-2*t)*(1-2*t+2*t**2))
def coeffs(e,N):
    s=sp.series(e,t,0,N).removeO(); p=sp.Poly(s,t)
    return [p.coeff_monomial(t**k) for k in range(N)]
def L(k,N):
    c=coeffs(h**k,N); return sp.Matrix(N,N,lambda i,r: c[i-r] if i>=r else 0)
def R(k,N): return sp.Matrix(N,N,lambda r,j: (-1)**(j-r)*sp.binomial(k,j-r) if 0<=j-r<=k else 0)
P=(3,0,3,3);B=(0,1,-2,-4);C=(1,3,-2,-8);D=(5,4,3,-5);E=(6,6,3,-9);Fm=(9,5,8,-2);H=(1,1,0,-2);I=(2,3,0,-6);J=(6,4,5,-3);K=(7,6,5,-7);T=(8,3,8,2);V=(5,2,5,1)
allok=True
def same_size(X,Y,k,m,detA0,detB0,Nmin,Nmax=11):
    global allok
    ok=True
    for N in range(Nmin,Nmax):
        M=L(k,N)*A(X,N)*R(k,N); AY=A(Y,N)
        ok&=all(M[i,j]==AY[i,j] for i in range(N) for j in range(m,N))
        ok&=all(AY[i,j]==0 for i in range(min(m,N)) for j in range(m,N))
        if N>=m: ok&=(M[:m,:m].det()==detA0 and AY[:m,:m].det()==detB0)
    print(X,'->',Y,'k',k,'m',m,'blocks',ok); allok&=ok
def shift(X,Y,c,Nmin,Nmax=11):
    global allok
    ok=True
    for N in range(Nmin,Nmax):
        M=L(1,N)*A(X,N)*R(1,N); AY=A(Y,N-1)
        ok&=M[0,0]==c and all(M[0,j]==0 for j in range(1,N)) and all(M[i,j]==AY[i-1,j-1] for i in range(1,N) for j in range(1,N))
    print(X,'=> shift',Y,'top',c,ok); allok&=ok
shift(P,T,2,2)
shift((2,1,2,0),(7,4,7,-1),2,2)
same_size(B,C,1,3,24,-10752,3)
same_size(P,D,2,2,16,1008,2)
same_size(P,E,3,3,1024,256512256,3)
same_size(T,Fm,1,1,8,-10,1)
same_size(H,I,1,2,-8,392,2)
same_size(V,J,1,1,4,28,1)
same_size(V,K,2,2,128,356928,2)
# Prop 22 pieces used
ok=all(G(B,n+1)==3*G(T,n-1) for n in range(1,10)) and all(G(H,n)==-2*G(V,n-1) for n in range(1,10)); print('Prop22 pieces',ok); allok&=ok
# final identities n=3..10
f=True
for n in range(3,11):
    v=[G(P,n),F(2,3)*G(B,n+1),F(-1,672)*G(C,n+1),F(1,63)*G(D,n),F(4,1002001)*G(E,n),F(-8,5)*G(Fm,n-1)]
    w=[G(H,n),F(-1,49)*G(I,n),F(-2,7)*G(J,n-1),F(-4,5577)*G(K,n-1)]
    f&=all(x==v[0] for x in v) and all(x==w[0] for x in w) and G((2,1,2,0),n)==2*G((7,4,7,-1),n-1)
print('identities (10.10)-(10.12) n=3..10',f); allok&=f
# uniform column identity, symbolic: regular u_j for (a,b,g,d)
X,Y,Z=sp.symbols('X Y Z',positive=True)  # X=(1-2t)^(-4j), Y=(1-t)^(-4j), Z=t^(4j)
def u(a,b,g,d,s):  # regular column gf at column j+s
    return 2**(b+a-g)*t**(a-g)*X*(1-2*t)**(-4*s-a-1)+(-1)**a*Z*t**(4*s+d+1)*Y*(1-t)**(-4*s-a-1)
uni=True
for (a,b,g,d) in [(3,0,3,3),(0,1,-2,-4),(1,1,0,-2),(5,2,5,1),(8,3,8,2),(2,1,2,0)]:
    lhs=h*(u(a,b,g,d,1)-u(a,b,g,d,0)); rhs=u(a+1,b+2,g,d-4,1)
    uni&=sp.simplify(lhs-rhs)==0
print('uniform difference identity (symbolic in j)',uni); allok&=uni
print('RESULT','PASS' if allok else 'FAIL')
