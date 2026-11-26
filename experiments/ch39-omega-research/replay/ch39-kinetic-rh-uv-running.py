"""Continue the threshold model above the crossing; bound its first one-loop inverse-coupling zero."""
import json
from pathlib import Path
from runpy import run_path

from flint import arb,arb_poly,ctx

directory=Path(__file__).parent
shell=run_path(str(directory/"ch39-kinetic-rh-gaussian-crossing.py"))


def finite_model(order):
    data=json.loads((directory/"ch39-kinetic-precision-data.json").read_text())
    row=data["models"][order-1]
    P,D=arb_poly(row["P"]),arb_poly(row["Q"])
    den=arb(2)**data["root_bits"]
    roots=[(arb(k)/den).union(arb(k+1)/den) for k in reversed(row["roots"])]
    return [-1/(1681*x) for x in roots],[-P(x)/(x*D.derivative()(x)) for x in roots]


def main(shared=False):
    ctx.prec,ctx.cap=3072,3
    base=shell["base"]
    a,b,c=[base["gram"](i,j) for i,j in ((1,0),(1,2),(0,2))]
    p=(1+2*a*b*c-a*a-b*b-c*c)/(1-c*c)
    modes=[(arb(6),7),(arb(4),6)]
    for j,n in ((5,3),(3,2)):
        modes += [(3-arb(j).sqrt(),n),(3+arb(j).sqrt(),n)]
    d=[8*arb.pi()**2/arb(v)**2 for v in (".467453",".63811","1.0583")]
    K=5*d[0]-3*d[1]-2*d[2]
    F=arb(109)/15*(d[1]-d[2])-arb(23)/6*(d[0]-d[1])
    Q=1000*(K/44).exp()
    cases=[(2,"3.5515295e13","3.5515297e13"),(21,"7.4768884e13","7.4768887e13"),(14,None,None)]
    if shared:
        cases=[(2,"2.4997146e13","2.4997148e13"),(14,"5.1757167e13","5.1757170e13"),
               (21,"5.9105125e13","5.9105128e13")]
    for order,lo,hi in cases:
        if order==14:
            z,w=finite_model(order)
        else:
            z,w,_=shell["quad"](order)
        mean=sum(z)/order
        C=[mean/v for v in z]
        rho=[(1-p)*a for a in w]
        assert min(rho)>0 and sum(w).contains(1) and all(u<v for u,v in zip(C,C[1:]))
        def masses(mf):
            spectrum=shell["scalar_spectrum"](mf,C,rho,p,modes)
            if shared:
                # The first m entries are the 25 cycle-only copies per pole.
                assert spectrum[:order]==[(v,25) for v in C]
                spectrum=spectrum[order:]
            return [(mf*v.sqrt(),n) for v,n in spectrum]
        def positive(t):
            if t>0:return t
            assert t<0
            return arb(0)
        def running(energy,mf,spec):
            E=sum(n*positive((energy/m).log()) for m,n in spec)
            trip=positive((energy/(mf*arb(8).sqrt())).log())
            sing=positive((energy/(mf*arb(24).sqrt())).log())
            T=(energy/1000).log()
            return [d[0]-arb(41)/10*T-E/10-arb(16)/5*trip-arb(4)/5*sing,
                    d[1]+arb(19)/6*T-E/6,d[2]+7*T-2*trip]
        if lo is None:
            # Bracket by integer-GeV bisection on the closed-fermion branch.
            left,right=20*10**12,200*10**12
            for _ in range(26):
                mid=(left+right)//2
                mf=arb(mid)
                E=sum(n*positive((Q/m).log()) for m,n in masses(mf))
                if 22*E>15*F:left=mid
                else:
                    assert 22*E<15*F
                    right=mid
            lo,hi=str(left),str(right)
        intervals=[]
        for sign,raw in [(1,lo),(-1,hi)]:
            mf=arb(raw)
            spec=masses(mf)
            E=sum(n*positive((Q/m).log()) for m,n in spec)
            assert sign*(22*E-15*F)>0 and Q/mf<arb(8).sqrt()
            active=sum(n for m,n in spec if m<Q)
            assert active>19 and min(m for m,n in spec)>1000
            poles=[]
            for group in (0,1):
                # Use a fixed finite positive energy bracket. At and above Q the
                # electroweak inverse slopes are negative; the color slope is positive.
                left,right=5*10**13,(10**27 if shared and order==2 else 10**21)
                assert running(arb(left),mf,spec)[group]>0 and running(arb(right),mf,spec)[group]<0
                for _ in range(72):
                    mid=(left+right)//2
                    val=running(arb(mid),mf,spec)[group]
                    if val.contains(0):
                        break
                    if val>0:left=mid
                    else:right=mid
                poles.append(arb(left).union(arb(right)))
            intervals.append(poles)
        U,V=[intervals[0][j].union(intervals[1][j]) for j in (0,1)]
        print(order,"reference modes:",active,"active at Q; M_F",lo,hi)
        print("One-loop inverse zeros for U(1), SU(2) in GeV",U.str(7),V.str(7))
        if shared:
            expected={2:(46,0,"6.5242e20","6.5245e20"),14:(95,1,"9.1848e15","9.1851e15"),
                      21:(97,1,"2.3494e15","2.3498e15")}
            count,group,low,high=expected[order]
            assert active==count and arb(low)<[U,V][group] and [U,V][group]<arb(high)
            continue
        if order==2:
            assert U>arb("1.3007e19") and U<arb("1.3009e19") and V>U
        if order==14:
            assert active==119 and V>arb("1.0507e15") and V<arb("1.0509e15") and U>V
        if order==21:
            assert V>arb("5.5348e14") and V<arb("5.5350e14") and U>V


if __name__=="__main__":
    main()
