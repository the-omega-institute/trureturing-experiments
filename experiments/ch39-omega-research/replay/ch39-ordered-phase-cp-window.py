"""The entire positive density range for the observed CP magnitude on the ordered branch."""
from pathlib import Path
import runpy
from flint import arb,acb_mat,arb_poly

phase=runpy.run_path(str(Path(__file__).with_name("ch39-symmetric-source-transition.py")))
fit=phase["fit"]
finite=phase["finite"]
_,wall,coeff=phase["data"](3)
B,s=phase["B"],phase["sym"]
I=finite["identity"]
adj=lambda m:m.conjugate().transpose()


def algebra(t):
    # Square the circulant background with the symmetric strength first, then insert
    # a rank-one change in each of three translated entries. Use the same mass matrix
    # up to a cyclic permutation in both sectors to avoid any eigenvector phase choices.
    k=acb_mat(wall)
    for i,w in enumerate((B-2*s*t,s*t,s*t)):
        finite["add_edge"](k,fit["inner"][i],fit["outer"][(i+1)%3],w)
    tr=-fit["block"](k,fit["ii"],fit["oo"])
    y=tr*tr
    return y*adj(y)


def trace(m):
    return sum(m[i,i] for i in range(3))


def invariants(t):
    h=algebra(t)
    d=acb_mat([[h[(i+1)%3,(j+1)%3] for j in range(3)] for i in range(3)])
    comm=h*d-d*h
    # Translate and scale the characteristic polynomial to avoid subtracting
    # large eigenvalue products. The discriminant is exactly the gap-product squared.
    scale=trace(h).real/3
    m=h-scale*I
    p=-trace(m*m).real/2
    q=-m.det().real
    disc=-4*p**3-27*q**2
    return comm.det().imag,2*disc


def certify():
    nodes=[arb(i)/24 for i in range(25)]
    v=[invariants(t) for t in nodes]
    p,d=[arb_poly.interpolate(nodes,values) for values in zip(*v)]
    # H is quartic in t and det[H,PHP*] and the discriminant are at most degree 24.
    deriv=p.derivative()*d-p*d.derivative()
    def verify(a,b):
        box=a.union(b)
        if p(box)>0 and d(box)>0 and deriv(box)<0:
            return 1
        assert b-a>arb("1e-7")
        mid=(a+b)/2
        return verify(a,mid)+verify(mid,b)
    cells=verify(arb(0),arb(".668"))
    a,c,f=coeff
    def tau(t):
        x=s*t
        y=B-2*x
        return x*y*(c+f*x)/(a+c*(x*x+2*x*y)+f*x*x*y)
    bounds=[]
    for target in (".0000325",".00003"):
        delta=p-arb(target)*d
        lo,hi=arb(0),arb(".667")
        assert delta(lo)>0>delta(hi)
        for _ in range(49):
            mid=(lo+hi)/2
            numerator,denominator=invariants(mid)
            value=numerator-arb(target)*denominator
            if value>0:
                lo=mid
            else:
                assert value<0
                hi=mid
        t=lo.union(hi)
        bounds.append(tau(t))
        print("PDG endpoint |J|, t, density exponent",target,t.str(13),tau(t).str(13))
    assert bounds[0]<bounds[1]<arb(".0000334")
    print("monotonicity boxes",cells)
    return bounds


if __name__=="__main__":
    certify()
