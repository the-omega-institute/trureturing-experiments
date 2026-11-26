"""Physical clock-unit bounds without normalizing the mean reference inertia."""
from pathlib import Path
from runpy import run_path

from flint import arb,ctx

directory=Path(__file__).parent
shell=run_path(str(directory/"ch39-kinetic-rh-gaussian-crossing.py"))


def main():
    ctx.prec,ctx.cap=1536,3
    z,w,_=shell["quad"](21)
    poles=[1/v for v in z]
    a,b,c=[shell["base"]["gram"](i,j) for i,j in ((1,0),(1,2),(0,2))]
    p=(1+2*a*b*c-a*a-b*b-c*c)/(1-c*c)
    rho=[(1-p)*v for v in w]
    modes=[(arb(6),7),(arb(4),6)]
    for j,n in ((5,3),(3,2)):
        modes.extend([(3-arb(j).sqrt(),n),(3+arb(j).sqrt(),n)])
    d=[8*arb.pi()**2/arb(v)**2 for v in (".467453",".63811","1.0583")]
    K=5*d[0]-3*d[1]-2*d[2]
    F=arb(109)/15*(d[1]-d[2])-arb(23)/6*(d[0]-d[1])
    A,B=arb(8).log()/2,arb(24).log()/2
    Q=1000*(K/44).exp()
    positive=lambda v:v if v>0 else arb(0) if v<0 else arb(0).union(abs(v))
    # Universal interlacing bound. The 23 lightest roots separately
    # dominate the static exposures; that branch gives the sharper w cap.
    m_static=arb("2034928.14").union(arb("2034928.16"))*(1/p).sqrt()
    q_static=arb("3.8981051e12").union(arb("3.8981054e12"))
    wcap=(q_static/m_static).log()
    required=(15*F+402*(wcap-A)+46*(wcap-B))/22
    def reference_floor(mu):
        return 48*sum(positive((q_static*z_i.sqrt()/mu).log()) for z_i in z)
    assert reference_floor(arb("1e10"))>required
    lo,hi=10**10,10**12
    for _ in range(50):
        if hi-lo<=1:break
        mid=(lo+hi)//2
        v=reference_floor(arb(mid))-required
        if v.contains(0):break
        if v>0:lo=mid
        else:hi=mid
    assert reference_floor(arb(lo))>required
    print("All physically allowed central-input crossings require clock unit >",lo,"GeV by interlacing.")
    # Finite reference roots at M_F growing: rank-one downward shift,
    # independently of lambda. The remaining one-per-block state decouples.
    limiting=[]
    for left,right in zip([arb(0),*poles[:-1]],poles):
        lo,hi=left.lower(),right.upper()
        for _ in range(92):
            mid=(lo+hi)/2
            val=p+mid*sum(r/(mid-u) for r,u in zip(rho,poles))
            if val.contains(0):break
            if val>0:lo=mid
            else:hi=mid
        root=lo.union(hi)
        assert left<root and root<right
        limiting.append(root)
    def load(masses):
        return sum(n*positive((Q/m).log()) for m,n in masses)
    for mu,sgn in [(arb("1.9612013e12"),1),(arb("1.9612015e12"),-1)]:
        spectrum=[(mu*v.sqrt(),25) for v in poles]+[(mu*v.sqrt(),23) for v in limiting]
        assert sgn*(22*load(spectrum)-15*F)>0
    mf=Q/arb(8).sqrt()
    for mu,sgn in [(arb("2.4033933e12"),1),(arb("2.4033937e12"),-1)]:
        spectrum=shell["scalar_spectrum"](mf,[(mu/mf)**2*v for v in poles],rho,p,modes)
        masses=[(mf*v.sqrt(),n) for v,n in spectrum]
        assert min(m for m,n in masses)>1000 and sgn*(22*load(masses)-15*F)>0
    # Lower corner for joint monotonicity of the physical mass branches.
    start=shell["scalar_spectrum"](mf,[(arb("1.9612013e12")/mf)**2*v for v in poles],rho,p,modes)
    assert min(mf*v.sqrt() for v,n in start)>1000
    print("Closed-fermion plateau clock-unit interval: lower endpoint in (1.9612013,1.9612015)e12, upper in (2.4033933,2.4033937)e12 GeV.")


if __name__=="__main__":
    main()
