"""One-loop crossing ranges over the full half-last-digit input box for the same shell record."""
from pathlib import Path
from runpy import run_path

from flint import arb, ctx

directory=Path(__file__).parent
older=run_path(str(directory/"ch39-finite-rh-propagating-references.py"))


def main():
    ctx.prec,ctx.cap=320,3
    gs=list(map(arb,(".467453",".63811","1.0583")))
    radii=list(map(arb,(".0000005",".000005",".00005")))
    dlow=[8*arb.pi()**2/(g+e)**2 for g,e in zip(gs,radii)]
    dhigh=[8*arb.pi()**2/(g-e)**2 for g,e in zip(gs,radii)]
    box=[v.union(w) for v,w in zip(dlow,dhigh)]
    keys=[5,-3,-2]
    fs=[-arb(23)/6,arb(109)/15+arb(23)/6,-arb(109)/15]
    dot=lambda a,b:sum(v*w for v,w in zip(a,b))
    a,b=arb(8).log()/2,arb(24).log()/2
    eig=[(arb(6),7),(arb(4),6)]
    for k,n in [(5,3),(3,2)]:
        eig += [(3-arb(k).sqrt(),n),(3+arb(k).sqrt(),n)]
    # All filtered or integer readings certified in the same six bins.
    R=arb("528.89869680").union(arb("528.89869682"))
    def scales(m,d):
        K,F=dot(keys,d),dot(fs,d)
        Z=sum(n*(v+125**2/m**2).log() for v,n in eig)/2
        A=K-44*(m/1000).log()+12*a+4*b
        B=15*F+22*Z-402*a-46*b
        ell=(448*A+16*B)/27808
        w=(506*A-44*B)/27808
        assert w>b and ell>2 and m**2*min(v for v,n in eig)+125**2>10**6
        return (ell-w).exp(),m*ell.exp()
    mlo,mhi=arb(1800000),arb(2300000)
    assert scales(mlo,box)[0]<R and scales(mhi,box)[0]>R
    mrange=mlo.union(mhi)
    H=sum(n*125**2/(125**2+mrange**2*v) for v,n in eig)
    rm=(2552-1320*H)/27808
    qm=1-(19712+352*H)/27808
    md=[(58*k-900*f)/(27808*rm) for k,f in zip(keys,fs)]
    qd=[(448*k+240*f)/27808+qm*m for k,f,m in zip(keys,fs,md)]
    assert rm>0 and qm>0 and all(v*w>0 for v,w in zip(md,qd))
    signs=[int(v.sgn().unique_fmpz()) for v in qd]
    assert signs==[1,-1,1]
    dl=[u if s>0 else v for s,u,v in zip(signs,dlow,dhigh)]
    dh=[v if s>0 else u for s,u,v in zip(signs,dlow,dhigh)]
    def bracket(d):
        low,high=1800000,2300000
        while high-low>1:
            mid=(low+high)//2
            r,_=scales(arb(mid),d)
            if r<R:
                low=mid
            else:
                assert r>R
                high=mid
        return arb(low).union(arb(high))
    lo,hi=bracket(dl),bracket(dh)
    qlo,qhi=scales(lo,dl)[1],scales(hi,dh)[1]
    assert qlo>arb("3.8670e12") and qhi<arb("3.9295e12")
    print("Static min and max M_S",lo.str(12),hi.str(12),"GeV; Q",qlo.str(10),qhi.str(10))

    # For the independent kinetic references the branch with both fermions
    # closed requires an exposure root at Q=1000 exp(K/44).
    q_min=1000*(dot(keys,[dlow[0],dhigh[1],dhigh[2]])/44).exp()
    q_max=1000*(dot(keys,[dhigh[0],dlow[1],dlow[2]])/44).exp()
    q=q_min.union(q_max)
    f=dot(fs,box)
    x,y,z=older["gram"](1,0),older["gram"](1,2),older["gram"](0,2)
    residual=(1+2*x*y*z-x*x-y*y-z*z)/(1-z*z)
    def exposure(mf):
        masses=[(1+z,25),(1-z,25)]
        for lam,num in eig:
            def poly(t):
                return (t-lam-(125/mf)**2)*((t-1)**2-z*z)-lam*((x*x+y*y)*(t-1)+2*x*y*z)
            for left,right in [("0",".0001"),(".001",".01"),("1","9")]:
                masses.append((older["root_bracket"](poly,left,right),num))
        total=arb(0)
        active=0
        for val,num in masses:
            assert mf*val.sqrt()>1000
            e=(q/(mf*val.sqrt())).log()
            if e>0:
                total+=num*e
                active+=num
            else:
                assert e<0
        assert active==71 and q/mf<arb(8).sqrt()
        return total
    assert 22*exposure(arb("2.1e14"))>15*f and 22*exposure(arb("2.4e14"))<15*f
    w0=arb(8).log()/2
    against=50*w0-arb(25)/2*(1-z*z).log()+23*(w0-arb(2).log()/2)
    against-=arb(23)/2*(arb(1)/64+6*residual/8).log()
    assert 22*against>15*f
    assert q_min>arb("4.99879e13") and q_max<arb("5.00478e13") and q_min/qhi>arb("12.72")
    print("Propagating-reference min and max Q",q_min.str(13),q_max.str(13),"GeV.")
    print("Independent-box lower bound for the ratio of crossing energies",(q_min/qhi).str(9))


if __name__=="__main__":
    main()
