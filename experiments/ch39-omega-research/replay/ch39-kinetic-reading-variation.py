"""Separation of critical readings for two exactly calibrated weights on the same 49 heights."""
from pathlib import Path
from runpy import run_path
from flint import arb

root=Path(__file__).resolve().parent
fiber=run_path(str(root/"ch39-vote-weighted-variation.py"))
base,ball=fiber["base"],fiber["ball"]


def main():
    poles,ends,difference_over_y=fiber["fiber"]()
    p,modes,q02,_=base["physical_inputs"]()
    p,q02=ball(p),ball(q02)
    c=arb(".00082657012849023").union(arb(".00082657012849025"))
    H=125**2*poles[2]/q02
    eta=arb(1)/128
    def slope(y):
        ds=[sum(w*v/(v-y)**2 for w,v in zip(row,poles)) for row in ends]
        return ds[0].union(ds[1])
    ec=et=arb(0)
    for lam,n in zip(map(ball,modes),base["MODES"]):
        def secular(y):
            vals=[c*(y-H)-lam*p+lam*(1-p)*y*sum(w/(v-y) for w,v in zip(row,poles)) for row in ends]
            return vals[0].union(vals[1])
        for j in range(3):
            start=poles[j-1].upper() if j else arb(0)
            stop=poles[j].lower()
            l,r=start,stop
            for _ in range(85):
                mid=(l+r)/2
                if secular(mid)<0:l=mid
                else:r=mid
            lo=l
            l,r=start,stop
            for _ in range(85):
                mid=(l+r)/2
                if secular(mid)>0:r=mid
                else:l=mid
            hi=r
            assert secular(lo)<0 and secular(hi)>0
            y=lo.union(hi)
            assert max(start,H)<y and y<stop
            fy=c+lam*(1-p)*slope(y)
            ec+=n*(1-H/y)/fy
            et+=n*lam*(1-p)*difference_over_y(y)/fy
    cprime=-et/ec
    y=-eta/c
    tprime=(difference_over_y(y)/c)-slope(y)*cprime/(c*c)
    assert arb("2.62661e-41")<tprime and tprime<arb("2.62664e-41")
    print("d c_star / dt",cprime.str(14),"d T_star / dt",tprime.str(14))


if __name__=="__main__":
    main()
