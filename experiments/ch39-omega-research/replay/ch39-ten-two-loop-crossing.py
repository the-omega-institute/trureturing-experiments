"""在部分两圈退耦处方的整类允许函数上，单三角刚度仍可实现精确交汇。"""
from pathlib import Path
from runpy import run_path

from flint import acb,arb,ctx

correction = run_path(str(Path(__file__).with_name("ch39-ten-light-two-loop-residual.py")))
vacuum,soft = correction["vacuum"],correction["soft"]


def verify():
    ctx.prec=320
    h=acb.zeta_zero(11431).imag
    r=soft["geo"]["theta"](804)/soft["geo"]["theta"](12)
    op,g,v,s,w=vacuum["setup"](h)
    kmin,kmax=arb(17500),arb(500000)
    radius=arb.pi()*(arb("28.5")/(2*kmin)).sqrt()
    first=(arb(3).sqrt()*radius).cos()
    assert first*kmin > 40
    local=g/(first*kmin-40)
    assert (arb(3).sqrt()*local).cos()>arb(".999") and v/kmin < local
    Lmin,Lmax=((arb(mu)/arb("91.1876")).log() for mu in ("1e15","1e16"))
    def family(K):
        x,displacement=vacuum["light"](K,op,v,s,w)
        return correction["insertion_family"](x,displacement,r)[0]
    # P eliminates the one-loop time slope. J is independent of P on the
    # quotient by the common reading. Opposite face signs give a joint root.
    for K,sign in ((kmin,-1),(kmax,1)):
        read=family(K)
        margins=[]
        for i in range(8):
            lo=Lmin+(Lmax-Lmin)*i/8
            hi=Lmin+(Lmax-Lmin)*(i+1)/8
            P=read(lo,hi,64,arb(3))[3][0]
            margins.append((sign*P).lower())
        print("K face signed P lower",K.str(8),min(margins).str(9))
        assert min(margins) > (arb(".03318") if sign == -1 else arb(".2529"))
    both=[]
    scale_margins = {1: [], -1: []}
    for j in range(16):
        a=kmin*(kmax/kmin)**(arb(j)/16)
        b=kmin*(kmax/kmin)**(arb(j+1)/16)
        read=family(a.union(b))
        for L,sign in ((Lmin,1),(Lmax,-1)):
            J=read(L,L,64,arb(3))[3][1]
            assert sign*J>0,(j,J)
            scale_margins[sign].append((sign*J).lower())
        both.append(read(Lmin,Lmax,64,arb(3))[3][2])
    for sign, margins in scale_margins.items():
        print("Scale face signed J lower", sign, min(margins).str(9))
        assert min(margins) > (arb("5.346") if sign == 1 else arb("2.583"))
    # A zero-sum projection costs at most its sum of positive weights times
    # the spread of a correction. For a hypercharge-only correction use
    # its individual P and J coefficients.
    assert arb(333)/218 * arb(".0217") < arb(".03318")
    assert 5 * arb(".0217") < arb("2.583")
    assert arb(115)/218 * arb(".0628") < arb(".03318")
    assert 5 * arb(".0628") < arb("2.583")
    common=min(z.lower() for z in both).union(max(z.upper() for z in both))
    print("Common inverse across the face-certified rectangle",common.lower().str(8),common.upper().str(8))
    assert arb("52.0") < common < arb("56.1")


if __name__=="__main__":
    verify()
