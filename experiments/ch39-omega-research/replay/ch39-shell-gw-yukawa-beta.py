"""A basis-invariant one-loop test of the fixed shell Yukawa tensor."""
from pathlib import Path
from runpy import run_path
from flint import acb_mat, arb

common=run_path(str(Path(__file__).with_name("ch39-shell-gw-195-common-field.py")))
adj,hs=common["adj"],common["hs"]


def blocks(F):
    return [acb_mat(3,3,[F[i,j] for i in range(9)]) for j in range(54)]


def trace(A):
    return sum(A[i,i] for i in range(A.nrows()))


def shape(F,other):
    y,z=blocks(F),blocks(other)
    G,H=adj(F)*F,adj(other)*other
    left=sum((a*adj(a)+b*adj(b) for a,b in zip(y,z)),acb_mat(3,3))/2
    right=sum((adj(a)*a for a in y),acb_mat(3,3))
    beta=[]
    for k in range(54):
        value=left*y[k]+y[k]*right
        for j in range(54):
            value+=3*(G[j,k]+H[k,j])*y[j]-2*z[j]*adj(z[k])*y[j]
        beta.append(value.entries())
    B=acb_mat(beta).transpose()
    A=F*adj(F)
    E=B*adj(F)+F*adj(B)
    trace0,trace2=trace(A).real,trace(A*A).real
    alpha=trace(E).real/(2*trace0)
    derivative=2*trace(A*E).real/trace0**2-2*trace2*trace(E).real/trace0**3
    return alpha,derivative


if __name__=="__main__":
    for n in (43,195):
        F=common["old"]["F"] if n==43 else common["operators"]()[1]
        t=common["native"]["relative_phase"](arb(676))
        G=common["conjugate"](F,t).conjugate()
        for a,b in ((F,G),(G,F)):
            alpha,slope=shape(a,b)
            assert arb("1.8745")<alpha<arb("1.8751")
            assert arb(".0343")<slope<arb(".0346")
        print("zero",n,"Yukawa cubic radial slope",alpha.str(10),"normalized shape drift",slope.str(10))
