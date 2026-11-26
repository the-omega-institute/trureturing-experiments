"""Sharp word-frame gain at prescribed low shortest-arc loss, checked with polynomial identities."""
from itertools import product
import sympy as sp

c,u,s,x=sp.symbols("c u s x",real=True)


def dictionary(n,phase):
    letters=[]
    for k in range(n):
        i,j=(k+1)%3,(k+2)%3
        v=sp.zeros(3);v[k,k]=1
        v[i,i],v[j,j]=c+sp.I*phase,c-sp.I*phase
        v[i,j],v[j,i]=-s,s
        letters.extend((v,v.adjoint()))
    return [sp.eye(3)]+letters+[v*w for v,w in product(letters,repeat=2)]


def normalized_polynomial(a):
    return sp.factor(sp.rem(sp.expand(a),s**2+c**2-1,s))


def nonnegative(a):
    expr=sp.cancel(a.subs(c,-(2+x)/(2+2*x)))
    top,bottom=sp.fraction(expr)
    assert all(v>=0 for v in sp.Poly(top,x).all_coeffs())
    assert all(v>=0 for v in sp.Poly(bottom,x).all_coeffs())
    assert a.subs(c,-sp.Rational(1,2))>=0


if __name__=="__main__":
    full=dictionary(3,u)
    rayleigh=sum((v[0,1]-v[1,0])*(v[0,1]-v[1,0]).conjugate() for v in full)/2
    rayleigh=sp.factor(sp.rem(sp.expand(rayleigh),s**2+c**2+u**2-1,s))
    assert sp.simplify(rayleigh-4*(1-c*c-u*u)*(5*c*c+4*c+4+u*u))==0
    for n in (2,3):
        words=dictionary(n,0)
        F=sp.Matrix([[w[i,j] for i,j in product(range(3),repeat=2)] for w in words])
        gram=(F.T*F).applyfunc(normalized_polynomial)
        target=4*(1-c*c)**2 if n==2 else 4*(1-c*c)*(5*c*c+4*c+4)
        groups=[[0,4,8],[1,3],[2,6],[5,7]]
        for j,g in enumerate(groups):
            for h in groups[j+1:]:
                assert gram.extract(g,h)==sp.zeros(len(g),len(h))
        diag=gram.extract(groups[0],groups[0])-target*sp.eye(3)
        for j in range(1,4):
            nonnegative(sp.factor(diag[:j,:j].det()))
        attained=False
        for i,j in groups[1:]:
            assert gram[i,i]==gram[j,j]
            for signed in (1,-1):
                gap=sp.factor(gram[i,i]+signed*gram[i,j]-target)
                nonnegative(gap)
                attained |= (gap==0)
        assert attained
    print("for -1 <= cos(p)cos(q) <= -1/2, exact upper gains attained at sin(p)=0")
