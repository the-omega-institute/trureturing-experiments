"""引入第二个物理弱双重态后，在三枚三代候选上重新反演强耦合。"""
from pathlib import Path
from runpy import run_path

from flint import acb, arb, ctx

ROOT = Path(__file__).parent
fixed = run_path(str(ROOT / "ch39-chiral-clock-unit-readout.py"))
chiral = fixed["chiral"]


def verify():
    model,_,_,initial,_ = chiral["source"]["mean"]["inputs"]()
    ctx.prec = 192
    read = chiral["source"]["spectral_supplier"](model)
    for n in (11661,11934,12207):
        rows = read(acb.zeta_zero(n).imag)
        layers = {q:chiral["ordered"]([h for h,z in rows if z==q])
                  for q in (0,*model["charges"])}
        assert all(layers[q][2]<9<layers[q][3] for q in model["charges"])
        print("零点与额外双重态质量/GeV",n,125.25)
        _,_,alpha = fixed["crossing_enclosure"](model,initial,layers,False,724,"28.947",".1252",
                                                additional_mass=arb("125.25"))
        assert arb(".1249")<alpha<arb(".1255")


if __name__=="__main__":
    verify()
