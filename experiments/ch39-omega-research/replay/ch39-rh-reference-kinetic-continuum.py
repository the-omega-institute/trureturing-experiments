"""Trace the crossing as both edge-reference kinetic coefficients grow from zero to one."""
from pathlib import Path
from runpy import run_path

from flint import acb_poly, arb, ctx

prior=run_path(str(Path(__file__).with_name("ch39-finite-rh-propagating-references.py")))


def main():
    ctx.prec=256
    d=[8*arb.pi()**2/arb(g)**2 for g in (".467453",".63811","1.0583")]
    K=5*d[0]-3*d[1]-2*d[2]
    F=arb(109)/15*(d[1]-d[2])-arb(23)/6*(d[0]-d[1])
    a,b,c=[prior["gram"](i,j) for i,j in ((1,0),(1,2),(0,2))]
    det=1+2*a*b*c-a*a-b*b-c*c
    pivot=det/(1-c*c)
    modes=[(arb(6),7),(arb(4),6)]
    for k,n in ((5,3),(3,2)):
        modes.extend([(3-arb(k).sqrt(),n),(3+arb(k).sqrt(),n)])
    A,B=arb(8).log()/2,arb(24).log()/2
    # At the static endpoint M_F>10^9; smaller w increases M_F.
    assert arb("2034928.14")*arb("528.89869680")>10**9
    for lam,n in modes:
        assert lam*det/(lam+2)**2*10**18>10**6
    assert (1-c)*10**18>10**6
    assert arb(1)/64+6*pivot/8<1
    assert 22*23*(1-arb(15)/(11*64))>448

    def boundary(w,tau,triplet,singlet):
        q=1000*((K-12*triplet-4*singlet)/44).exp()
        mf=q/ w.exp()
        z=(125/mf)**2
        minus,plus=(1-c)/tau,(1+c)/tau
        spectrum=[(minus,25),(plus,25)]
        X=acb_poly([0,1])
        for lam,n in modes:
            vm=lam*(a-b)**2/(2*tau)
            vp=lam*(a+b)**2/(2*tau)
            poly=(X-lam-z)*(X-minus)*(X-plus)-vm*(X-plus)-vp*(X-minus)
            roots=sorted((r.real for r in poly.roots(1e-14)),key=lambda t: float(t.mid()))
            l,m,h=roots
            assert 0<l and l<minus and minus<m and m<plus and plus<h
            spectrum.extend((val,n) for val in roots)
        exposure=arb(0)
        active=0
        for val,n in spectrum:
            term=w-val.log()/2
            if term>0:
                exposure+=n*term
                active+=n
            else:
                assert term<0
        assert active==71
        return 22*exposure-402*triplet-46*singlet-15*F,q

    for w,trip,sing,lo,hi in [
        (B,B-A,arb(0),".00027313416",".00027313419"),
        (A,arb(0),arb(0),".00091244947",".00091244951")
    ]:
        low,q=boundary(w,arb(lo),trip,sing)
        high,_=boundary(w,arb(hi),trip,sing)
        assert low<0 and high>0
        print("Reference kinetic threshold bracket",lo,hi,"with Q",q.str(14),"GeV.")
    print("The crossing spans the whole interval between the static energy and the closed-fermion plateau.")


if __name__=="__main__":
    main()
