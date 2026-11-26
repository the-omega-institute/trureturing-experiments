"""First arithmetic transition after charging the squared inverse Wilson-frame gain."""
from itertools import product
from pathlib import Path
from runpy import run_path

import numpy as np
from flint import acb, acb_mat, arb,ctx

w=run_path(str(Path(__file__).with_name("ch39-shell-gw-wilson-coordinates.py")))
ctx.prec=160
f,g=[a-w["angles"][0] for a in w["angles"][1:]]
scale=arb(1)/51**5


def readings(height):
    p,q=f*height,g*height
    return 1+p.cos()*q.cos()+scale*height**2,p,q


def gain(p,q,r):
    gens=[]
    for k in range(r):
        i,j=(k+1)%3,(k+2)%3
        U=acb_mat(3,3)
        U[k,k]=1
        U[i,i]=acb(p.cos(),p.sin())*q.cos()
        U[j,j]=U[i,i].conjugate()
        U[i,j],U[j,i]=-q.sin(),q.sin()
        gens.extend((U,w["adj"](U)))
    form=w["frame"]([w["identity"]]+gens+[a*b for a,b in product(gens,repeat=2)])
    rows=[[form[i,j] for j in range(9)] for i in range(9)]
    numeric=float(np.linalg.eigvalsh(np.array([[complex(z) for z in row] for row in rows]))[0])
    radius=arb("2e-13")
    lo,hi=arb(numeric)-radius,arb(numeric)+radius
    assert lo>0 and w["count"](rows,lo)==0 and w["count"](rows,hi)>0
    return lo.union(hi)


def switch(r,a,b,zeros):
    # a,b contain the uncharged loss and both phase coordinates.
    g0,g1=gain(*a[1:],r),gain(*b[1:],r)
    zeta=(b[0]-a[0])/(1/g0-1/g1)
    upper=zeta*arb("1.02")
    ceiling=a[0]+upper/g0
    assert scale*2300**2>ceiling
    for j,(loss,p,q) in enumerate(zeros,start=1):
        if j in (43,195):
            continue
        assert loss>a[0]
        if not loss>ceiling:
            val=gain(p,q,r)
            assert loss+upper/val>ceiling
    print("Wilson loops",r,"first penalty switch interval",zeta.str(13),
          "then zero 195 through",upper.lower().str(10))
    return zeta


if __name__=="__main__":
    total=int(arb(2300).zeta_nzeros().unique_fmpz())
    zeros=[]
    for start in range(1,total+1,128):
        zeros.extend(readings(z.imag) for z in acb.zeta_zeros(start,min(128,total+1-start)))
    assert len(zeros)==total
    z2=switch(2,zeros[42],zeros[194],zeros)
    z3=switch(3,zeros[42],zeros[194],zeros)
    assert arb(".00000058657027")<z2<arb(".00000058657031")
    assert arb(".00091367607")<z3<arb(".00091367610")
    print("all competing positive heights below 2300:",total)
    chiral=run_path(str(Path(__file__).with_name("ch39-shell-gw-selected-zero.py")))
    ctx.prec=288
    height=acb.zeta_zero(195).imag
    for sign in (1,-1,0):
        base,k=chiral["action_kernel"](sign,chiral["f"]*height,chiral["g"]*height)
        count=chiral["model"]["model"]["count"]
        gap=arb(".005")
        assert count(k,-gap)==count(k,gap)==18-3*sign
        chiral["model"]["calc"]["opposite_modes"](k,3*sign,gap)
        for sector in chiral["character_blocks"](base,k):
            assert count(sector,-gap)==count(sector,gap)==6-sign
    print("zero 195 retains three positive, three negative, and no neutral zero modes")
