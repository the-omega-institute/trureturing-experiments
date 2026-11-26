"""用一枚有编号的高处零点检验二十七圈作用的自由真空能量差。"""
from pathlib import Path
from runpy import run_path
from flint import acb, arb, ctx


def verify():
    root = Path(__file__).resolve().parent
    vacuum = run_path(str(root / "ch39-mincut-global-winding.py"))
    physical = run_path(str(root / "ch39-first-three-mode-observables.py"))
    shell = run_path(str(root / "ch39-three-ray-chiral-boundary.py"))
    ctx.prec = 144
    index = 4_224_080_297_279
    height = acb.zeta_zero(index).imag
    assert height > 10**12
    ctx.prec = 224
    phases = shell["geo"]["phases"](height)
    x, y = [phases[d].arg() for d in (724,804)]
    model, initial = physical["model_inputs"]()
    chiral = physical["joint"]["chiral"]
    rows = chiral["source"]["spectral_supplier"](model)(height)
    bands = {q:chiral["ordered"]([v for v,charge in rows if q==charge])
             for q in (0,*model["charges"])}
    delta = arb(23)/25*arb(2)**-20/(height**2+shell["soft"]["point"]["study"]["offset"])
    threshold = 9+(1-delta)**2
    assert all(bands[q][2] < threshold < bands[q][3] for q in model["charges"])
    _, alpha, curve = physical["layer_observables"](model, initial, bands)
    assert arb(".117")<alpha<arb(".119") and arb(".99953087404")<curve<arb(".99953087407")
    triangles = vacuum["terms"](shell)
    modes = [(q,w,vacuum["coeff"](triangles,q,delta)) for q,w in ((1,100),(27,13))]
    def action(a,b):
        return sum(w*vacuum["energy"](q,c,a,b) for q,w,c in modes)/113

    half = arb(0,"2.1e-17")
    vacuum_x, vacuum_y = [arb(c)+half for c in (".0707075643615173",".07570197643661009")]
    gap = action(x,y)-action(vacuum_x,vacuum_y)
    assert arb("2.2647298e-7") < gap < arb("2.2647300e-7")
    print("编号、高度、双钟相位", index, height.str(23), x.str(18), y.str(18))
    print("强耦合、四层曲率比、同高度真空能差", alpha.str(12), curve.str(15), gap.str(12))


if __name__ == "__main__":
    verify()
