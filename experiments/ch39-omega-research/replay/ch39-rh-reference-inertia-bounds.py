"""Sharp threshold envelopes with arbitrary common two-reference kinetic form."""
from pathlib import Path
from runpy import run_path

from flint import arb,ctx

base=run_path(str(Path(__file__).with_name("ch39-finite-rh-propagating-references.py")))


def main():
    ctx.prec=256
    d=[8*arb.pi()**2/arb(g)**2 for g in (".467453",".63811","1.0583")]
    K=5*d[0]-3*d[1]-2*d[2]
    F=arb(109)/15*(d[1]-d[2])-arb(23)/6*(d[0]-d[1])
    a,b,c=[base["gram"](i,j) for i,j in ((1,0),(1,2),(0,2))]
    p=(1+2*a*b*c-a*a-b*b-c*c)/(1-c*c)
    modes=[(arb(6),7),(arb(4),6)]
    for k,n in ((5,3),(3,2)):
        modes+= [(3-arb(k).sqrt(),n),(3+arb(k).sqrt(),n)]
    A,B=arb(8).log()/2,arb(24).log()/2
    def setup(w,trip,sing):
        return 1000*((K-12*trip-4*sing)/44).exp()/w.exp()
    def load(w,t,mf):
        e=arb(0)
        active=0
        def add(val,n):
            nonlocal active,e
            term=w-val.log()/2
            if term>0:
                e+=n*term
                active+=n
            else:
                assert term<0
        add(1/t,73)
        for lam,n in modes:
            z=(125/mf)**2
            tr=lam+z+1/t
            det=(lam*p+z)/t
            large=(tr+(tr*tr-4*det).sqrt())/2
            add(det/large,n)
            add(large,n)
        return e,active
    for w,trip,small,large,trans_lo,trans_hi,num in [
        (B,B-A,".05989035",".05989036",".08608450",".08608451",119),
        (A,arb(0),".17995156",".17995158",".28776332",".28776335",101)]:
        mf=setup(w,trip,0)
        target=(15*F+402*trip)/22
        lo,nlo=load(w,arb(small),mf)
        hi,nhi=load(w,arb(large),mf)
        assert nlo==nhi==num and lo<target and hi>target
        static=sum(n*(w-(lam*p+(125/mf)**2).log()/2) for lam,n in modes)
        assert all(w>(lam*p+(125/mf)**2).log()/2 for lam,_ in modes)
        transverse=(2*((target-static)/48-w)).exp()
        assert arb(trans_lo)<transverse and transverse<arb(trans_hi)
        print("At squared fermion threshold",w.exp()**2,"isotropic first reaches closure at",small,large,
              "all orientations reach closure by",transverse.str(11))
    q=1000*(K/44).exp()
    for mf,sign in [(arb("3.8101327e13"),1),(arb("3.8101328e13"),-1)]:
        w=(q/mf).log()
        exposure,n=load(w,arb(1),mf)
        assert sign*(22*exposure-15*F)>0 and n==96 and w<A
    print("Gram reference kinetic form at unit normalization: Q on plateau, 96 extra scalars open, M_F near 3.81013273e13 GeV.")


if __name__=="__main__":
    main()
