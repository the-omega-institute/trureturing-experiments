"""Joint minimization over the zero index and the shared continuous phase of the new shells."""
from pathlib import Path
from runpy import run_path
from flint import arb, acb

model = run_path(str(Path(__file__).with_name("ch39-shell-continuous-extension.py")))
pi = arb.pi()
penalty, choices = model["penalty"], model["choices"]
zero, point, refs = model["zero"], model["point"], model["refs"]


def dot(a,b):
    return sum(x*y for x,y in zip(a,b))


vectors={}
for r in (26,29):
    at0, at90, at180 = [refs(p)[r] for p in (0,arb(9)/4,arb(9)/2)]
    u=[(v+w)/2 for v,w in zip(at0,at180)]
    v=[x-y for x,y in zip(at0,u)]
    w=[x-y for x,y in zip(at90,u)]
    vectors[r]=u,v,w


def series(layers,hstate):
    base, measured = hstate
    coeff=[base,arb(0),arb(0),arb(0),arb(0)]
    for r in layers:
        u,v,w=vectors[r]
        z=[x-y for x,y in zip(measured[r],u)]
        block=[dot(z,z)+(dot(v,v)+dot(w,w))/2, -2*dot(z,v), -2*dot(z,w),
               (dot(v,v)-dot(w,w))/2,dot(v,w)]
        coeff=[x+y for x,y in zip(coeff,block)]
    return coeff


def value(coeff,x):
    # Here x is the reference angle divided by pi.
    s,c=(2*pi*x).sin_cos()
    t,b=(4*pi*x).sin_cos()
    return coeff[0]+coeff[1]*c+coeff[2]*s+coeff[3]*b+coeff[4]*t


def derivative(co):
    # Derivative with respect to 2*pi*x.
    return [arb(0),co[2],-co[1],2*co[4],-2*co[3]]


def minima(co):
    first=derivative(co)
    second=derivative(first)
    roots=[]
    def isolate(lo,hi):
        h=lo.union(hi)
        if not value(first,h).contains(0) or value(second,h)<0:
            return
        if value(second,h)>0:
            a,b=value(first,lo),value(first,hi)
            if a>0 or b<0:
                return
            assert a<0<b
            for _ in range(55):
                mid=(lo+hi)/2
                val=value(first,mid)
                if val<0:
                    lo=mid
                else:
                    assert val>0
                    hi=mid
            box=lo.union(hi)
            roots.append((value(co,box),box))
            return
        assert hi-lo>arb(1)/10000000
        mid=(lo+hi)/2
        isolate(lo,mid)
        isolate(mid,hi)
    for j in range(16):
        isolate(arb(j)/16,arb(j+1)/16)
    assert roots
    best=min(roots,key=lambda t:float(t[0].mid()))
    assert all(e>best[0] for e,x in roots if x is not best[1])
    return best


def main():
    candidates={rs:{719:minima(series(rs,zero(719)))} for rs in choices}
    upper={rs:row[719][0].upper() for rs,row in candidates.items()}
    cap=max(int((u/penalty).sqrt().upper().ceil().unique_fmpz()) for u in upper.values())
    assert all(penalty*cap**2>u for u in upper.values())
    # The second readout of each extra shell is independent of the freely moving first angle.
    ref=refs(0)
    def scan(lo,hi,parts):
        base,obs=point(lo.union(hi))
        keep=[rs for rs in parts if not base+sum((obs[r][1]-ref[r][1])**2 for r in rs)>upper[rs]]
        if not keep:
            return
        if hi-lo>arb(1)/4:
            mid=(lo+hi)/2
            scan(lo,mid,keep)
            scan(mid,hi,keep)
            return
        a,b=(v.zeta_nzeros().unique_fmpz() for v in (lo,hi))
        assert a is not None and b is not None
        for n in range(int(a)+1,int(b)+1):
            if n==719:
                continue
            for rs in keep:
                entry=minima(series(rs,zero(n)))
                candidates[rs][n]=entry
                if entry[0]<upper[rs]:
                    upper[rs]=entry[0].upper()
    for h in range(cap):
        scan(arb(h),arb(h+1),choices)
    for rs,row in candidates.items():
        n=min(row,key=lambda n:float(row[n][0].mid()))
        v,p=row[n]
        assert n==719
        assert all(z>v for k,(z,_) in row.items() if k!=n)
        print("shells",rs,"winner",n,"ref/pi",p.str(14),"min",v.str(14),"719",row[719][0].str(14))
    print("initial uniform tail",cap)


if __name__=="__main__":
    main()
