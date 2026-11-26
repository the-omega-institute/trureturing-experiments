"""Distinguish the two clocks by the intensity-weighted mean of the weak mass lines."""
from pathlib import Path
import runpy
from flint import acb_mat

m=runpy.run_path(str(Path(__file__).with_name("ch39-clock-weak-cubic.py")))
q,short=m["q"],m["short"]
u,v=m["assign"][7][1],m["assign"][8][1]

def phases(scales,derivatives):
    H=acb_mat(6,6)
    deriv=[acb_mat(6,6) for _ in range(4)]
    for _,(ell,g,D,A) in m["blocks"].items():
        y=ell*scales[g]
        sn,cs=y.sin_cos()
        b=y/sn
        H+=D*(b*cs)-A*b
        if derivatives:
            deriv[g]+=(D*((cs-b)/sn)-A*((1-b*cs)/sn))*ell
    small=acb_mat([[H[a,b] for b in range(1,6)] for a in range(1,6)])
    rest=small.solve(acb_mat([[-H[a,0]] for a in range(1,6)]))
    inner=rest[1,0]
    if inner.real.contains(0):
        return None
    slope=inner.imag/inner.real
    if not derivatives:
        return slope
    z=acb_mat([[1]]+[[rest[a,0]] for a in range(5)])
    rhs=[-d*z for d in deriv]
    dz=small.solve(acb_mat([[rhs[g][a,0] for g in range(4)] for a in range(1,6)]))
    return [(dz[1,g].imag-slope*dz[1,g].real)/inner.real for g in range(4)]

def ratio(box):
    f=short["lengths"](box[:3])
    if not f[3].overlaps(short["fsbox"][3]):
        return None
    f[3]=f[3].intersection(short["fsbox"][3])
    p=[box[3]*val/2 for val in f]
    mid=[x.mid() for x in p]
    gradient=phases(p,True)
    if gradient is None:
        return q(-2).union(q(0))
    phase=phases(mid,False)+sum((x-z)*d for x,z,d in zip(p,mid,gradient))
    # 724 and 804 have the same length multiplier. Their common scale
    # cancels in the prefactor of the ratio of endpoint weights.
    l,g,_,_=m["blocks"][724]
    h,j,_,_=m["blocks"][804]
    assert j==g
    x=p[g]
    b=(h/l)*(x*l).sin()/(x*h).sin()
    return b*(v.imag+v.real*phase)/(u.imag+u.real*phase)

def check():
    base=[*short["fsbox"][:3],q("1.19576").union(q("1.21744"))]
    stack=[(base,0)]
    values=[]
    while stack:
        box,depth=stack.pop()
        z=ratio(box)
        if z is None:
            continue
        if q("-1.10") < z < q("-1.010"):
            values.append(z)
            continue
        assert depth<30
        i=(1,2,3,0)[depth%4]
        for part in m["split"](box[i],2):
            sub=box[:]
            sub[i]=part
            stack.append((sub,depth+1))
    bound=min(x.lower() for x in values).union(max(x.upper() for x in values))
    delta=q(".00617")*(bound.upper()**2-1)
    at_target=delta*(m["s"]-m["t"])
    print("signed current ratio J804/J724",bound.lower().str(12),bound.upper().str(12),"boxes",len(values))
    print("minimum difference of oscillator shares",delta.lower().str(12))
    print("minimum 719 centroid change after transposing clocks",at_target.lower().str(12))

if __name__=="__main__":
    check()
