"""Certified switch of the fully relaxed phase model between the 255th and 25th zeros."""
from pathlib import Path
from runpy import run_path
from flint import arb,acb,arb_mat

home=Path(__file__).resolve().parent
m=run_path(str(home/"ch39-shell-four-free-references.py"))
pi=arb.pi()
idx=range(3)


def differentiate(obs,x):
    # Derivatives in radians; x is supplied in multiples of pi as in the parent model.
    H=arb_mat(3,3)
    grad=arb_mat(3,1)
    def add(res,variables,v,W):
        for i,k in enumerate(variables):
            grad[k,0]+=2*res*v[i]
            for j,l in enumerate(variables):
                H[k,l]+=2*(v[i]*v[j]+res*W[i][j])
    for r,(i,j) in m["layout"].items():
        c=m["m"]["shells"][r][1]
        s,b=(pi*x[i]).sin_cos()
        t,d=(pi*x[j]).sin_cos()
        f=b*d+c*s*t
        vp=-s*d+c*b*t
        vq=-b*t+c*s*d
        cross=s*t+c*b*d
        add(s*s-obs[r][0],[i],[2*s*b],[[2*(b*b-s*s)]])
        add(t*t-obs[r][1],[j],[2*t*d],[[2*(d*d-t*t)]])
        h12=-2*(vp*vq+f*cross)
        add(1-f*f-obs[r][2],[i,j],[-2*f*vp,-2*f*vq],
            [[2*(f*f-vp*vp),h12],[h12,2*(f*f-vq*vq)]])
    return grad,H


def pos(H):
    return H[0,0]>0 and H[0,0]*H[1,1]>H[0,1]*H[1,0] and H.det()>0


eye=arb_mat([[int(i==j) for j in idx] for i in idx])


def minimum(n,proposal):
    h=acb.zeta_zero(n).imag
    obs=m["state"](h)
    x=[arb(t) for t in proposal]
    for _ in range(3):
        g,H=differentiate(obs,x)
        correction=H.solve(g)
        x=[(t-correction[i,0]/pi).mid() for i,t in enumerate(x)]
    energy=m["losses"](obs,x)
    g,H=differentiate(obs,x)
    w=arb(1)/512
    local=[(v-w).union(v+w) for v in x]
    _,hess=differentiate(obs,local)
    assert pos(hess-eye/4)
    lower=energy-2*sum(g[i,0]**2 for i in idx)
    leaves=0
    def cover(box):
        nonlocal leaves
        if all(local[i].contains(lo.union(hi)) for i,(lo,hi) in enumerate(box)):
            leaves+=1
            return
        if m["interval_lower"](obs,box)>energy:
            leaves+=1
            return
        i=max(idx,key=lambda k:float(box[k][1]-box[k][0]))
        a,b=box[i]
        assert b-a>arb(1)/2**20
        mid=(a+b)/2
        for interval in ((a,mid),(mid,b)):
            cover(box[:i]+(interval,)+box[i+1:])
    cover(((arb(0),arb(1)/2),(arb(0),arb(1)),(arb(0),arb(1))))
    result=lower.union(energy)
    assert result.rad()<arb("1e-20")
    print("zero",n,"unique minimum modulo simultaneous sign, stationary angle estimates /pi",*[t.str(12) for t in x],
          "geometric loss",result.str(17),"boxes",leaves)
    return h,result


h25,e25=minimum(25,(".45246442",".56234423",".11710930"))
h255,e255=minimum(255,(".16209666",".54605929",".13268280"))
lam=arb(1)/51**5
cross=(e25-e255)/(h255*h255-h25*h25)
assert lam<cross<2*lam
# Any n>255 loses to 255 at lam by the earlier full-spectrum certificate and thereafter,
# since the coefficient of lambda is larger. At n<=255 test both ends against 25.
for n,z in enumerate(acb.zeta_zeros(1,255),1):
    if n in (25,255):
        continue
    h=z.imag
    low=m["paired_lower"](m["state"](h))
    for v in (lam,2*lam):
        assert low+v*h*h>e25+v*h25*h25
print("sharp switch inside [lambda,2 lambda]",cross.str(17),"multiple of nominal",(cross/lam).str(16))
