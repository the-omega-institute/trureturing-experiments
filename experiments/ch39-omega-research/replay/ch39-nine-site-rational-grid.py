"""十三个真实壳的弱差对九站有理门网格的统一有限误差。"""
from pathlib import Path
from runpy import run_path
from flint import arb, ctx


def verify(root, levels):
    def lib(name):
        return run_path(str(root/("ch39-"+name+".py")))
    ctx.prec=256
    floor, _, top=levels
    mz=arb("91.1876")
    rho=(arb("246.22")/mz)**2/2
    w=arb(".63811")**2*rho/2
    z=w+3*arb(".467453")**2*rho/10
    t=arb(".8616")**2*rho
    endpoint=(w/2+z/4+arb(1)/32+3*t/8+arb(1)/16)/arb.pi()
    graph=lib("three-ray-chiral-boundary")
    native=lib("native-rh-ray-pvm")
    gram=native["gram_cholesky"](root)[3]
    worst=arb(0)
    for j in range(7,20):
        _,basis=lib("adjacent-pvm-transport")["spectrum"](graph,6,2_306_888_100_085+j)
        weights=native["readout"](basis,range(19),[14,8,7,6],gram,3)
        variation=sum(abs(value.real).upper() for value in weights)
        worst=max(worst,(90*4*endpoint*((top-floor)/floor).sqrt()*variation).upper())
    for bits in (67,452):
        print("半角网格位数、十三壳最大弱差绝对误差上界",bits,(worst/2**bits).str(13),flush=True)
    assert worst/2**67 < arb("6e-7") and worst/2**452 < arb("7.6e-123")


if __name__=="__main__":
    root=Path(__file__).resolve().parent
    original=run_path(str(root/"ch39-fifth-full-positive-joint.py"))["verify"]
    verify(root,original(root,return_operator="spectral")[2])
