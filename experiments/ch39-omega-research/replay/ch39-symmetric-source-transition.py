"""A cyclic source with a permutation-invariant logarithmic density has an ordered phase."""
from pathlib import Path
import runpy
from flint import arb,acb,acb_mat,arb_poly

finite=runpy.run_path(str(Path(__file__).with_name("ch39-finite-orthogonal-ckm.py")))
fit=finite["fit"]
isolate=runpy.run_path(str(Path(__file__).with_name("ch39-orbit-flux-windows.py")))["isolate"]
B=arb(".1")
sym=B/3


def data(orbit):
    h=acb.zeta_zero(10).imag
    m,orbits=fit["model"](*(fit["s"][key]*h for key in ("a","b")))
    for v in orbits[orbit]:
        m+=arb(".05")*v*v.conjugate().transpose()
    wall,J,block=fit["W"],fit["J"],fit["block"]
    k=block(m,wall,wall)-block(m,wall,J)*block(m,J,J).solve(block(m,J,wall))
    H=acb_mat(19,3)
    for i in range(3):
        a,b=fit["inner"][i],fit["outer"][(i+1)%3]
        H[a,i]=1
        H[b,i]=-finite["bond"](a,b).conjugate()
    r=H.conjugate().transpose()*m.solve(H)
    coeff=1+B*r[0,0].real,r[0,0].real*r[1,1].real-abs(r[0,1])**2,r.det().real
    assert all(w>0 for w in coeff)
    return m.det().real,k,coeff


def classify():
    mats=[data(i) for i in range(4)]
    def score(coeff,x):
        a,c,d=coeff
        y=B-2*x
        return a+c*(x*x+2*x*y)+d*x*x*y
    det,k,coeff=mats[3]
    assert all(mats[i][0]*mats[i][2][0]>det*score(coeff,sym)*arb("1.00001") for i in range(3))
    a,c,d=coeff
    # In the dimensionless coordinate t=3x/B, all candidate minima have 0<t<=1.
    x=arb_poly([0,sym])
    y=B-2*x
    P=a+c*(x*x+2*x*y)+d*x*x*y
    N=x*y*(c+d*x)
    extrema=isolate(N.derivative()*P-N*P.derivative())
    assert len(extrema)==1 and 0<extrema[0]<1
    peak=extrema[0]
    def tau(t):
        return N(t)/P(t)
    def difference(t):
        return (P(t)/P(1)).log()-tau(t)*(t*t*(3-2*t)).log()
    # Along the stable unequal branch the derivative of the energy gap w.r.t. tau
    # is -log(t^2(3-2t))>0; comparing endpoint signs isolates coexistence.
    lo,hi=arb(".1"),arb(".72")
    assert hi<peak and difference(lo)<0<difference(hi)
    for _ in range(65):
        mid=(lo+hi)/2
        if difference(mid)>0:
            hi=mid
        else:
            assert difference(mid)<0
            lo=mid
    co=lo.union(hi)
    transition=tau(co)
    assert tau(1)<transition<tau(peak)
    def matrix(t):
        v=acb_mat(k)
        for i,w in enumerate((B-2*sym*t,sym*t,sym*t)):
            finite["add_edge"](v,fit["inner"][i],fit["outer"][(i+1)%3],w)
        return v
    # Off-cut verification on the entire one-large/two-equal segment.
    wall=fit["W"]
    def open_sheet(a,b):
        m=matrix(a.union(b))
        for face in fit["src"]["all_faces"]():
            z=acb(1)
            for i,j in zip(face,face[1:]+face[:1]):
                z*=-m[wall.index(i),wall.index(j)]
            if z.real>0 or not z.imag.contains(0):
                continue
            mid=(a+b)/2
            assert b-a>arb("1e-9")
            return open_sheet(a,mid)+open_sheet(mid,b)
        return 1
    cells=open_sheet(arb(0),arb(1))
    assert finite["flux"]["degrees"](matrix(0))[0]==(3,[3])
    roots=isolate(N-arb(".00001")*P)
    assert len(roots)==1 and roots[0]<co
    t=roots[0]
    def squared(mat):
        tr=-fit["block"](mat,fit["ii"],fit["oo"])
        return tr*tr
    vals=finite["observables"](squared(matrix(1)),squared(matrix(t)))
    print("two stability limits, coexistence, t at jump",
          tau(1).str(13),tau(peak).str(13),transition.str(13),co.str(13))
    print("flux interval cells, weak and strong stiffness at tau=1e-5, J",
          cells,(sym*t).str(13),(B-2*sym*t).str(13),vals[-1].str(12))
    assert arb("-.00003701")<vals[-1]<arb("-.00003700")
    def displaced(t,position):
        v=acb_mat(k)
        for i in range(3):
            strength=B-2*sym*t if i==position else sym*t
            finite["add_edge"](v,fit["inner"][i],fit["outer"][(i+1)%3],strength)
        return squared(v)
    js=[]
    for position in (1,2):
        val=finite["observables"](displaced(t,0),displaced(t,position))
        js.append(val[-1])
        assert all(arb(".044")<r<arb(".056") for r in val[:4])
    assert js[0]<arb("-.0000415") and js[1]>arb(".0000415")
    print("two ordered sectors, signed J for the two cyclic offsets",*(v.str(11) for v in js))
    return transition


if __name__=="__main__":
    classify()
