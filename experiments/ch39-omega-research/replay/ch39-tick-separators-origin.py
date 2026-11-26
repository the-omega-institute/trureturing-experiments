"""Build the finite inputs for the integer separation test from certified RH height cells."""
import json
from pathlib import Path
from runpy import run_path
from flint import acb, arb, arb_mat, ctx


def floor_at(x,bits):
    value=(x*2**bits).floor().unique_fmpz()
    assert value is not None
    return int(value)


def data():
    root=Path(__file__).resolve().parent
    five=json.loads((root/"ch39-native-five-tick-boundary.json").read_text())
    other=json.loads((root/"ch39-native-eight-tick.json").read_text())
    ctx.prec,ctx.cap=1792,3
    dual={}
    for tick in range(1,8):
        rec,depth=(five,48) if tick==5 else (other[str(tick)],44)
        z=[1/acb.zeta_zero(i).imag**2 if i else arb(0) for i in rec["nodes"]]
        d=len(z)-1
        v=arb_mat([[1/(1+j*j*t) for j in range(d+1)] for t in z])
        a=v.solve(arb_mat([[arb(n)/2**depth] for n in rec["values"][:-4]])).entries()
        dual[str(tick)]={"a":[floor_at(t,40) for t in a],"c":rec["values"][-4:],"depth":depth}
    ctx.prec=544
    heights=[v.imag for v in acb.zeta_zeros(1,3200)]
    logs=[arb(v).log() for v in (2,3,arb(3)/2,6)]
    return {"dual":dual,"prefix":[[floor_at(h,400),*[floor_at((h*x).cos(),72) for x in logs]] for h in heights]}


if __name__=="__main__":
    root=Path(__file__).resolve().parent
    assert data()==json.loads((root/"ch39-tick-separators.json").read_text())
    print("3,200 height cells, four base cosine cells per height, and seven rational separators certified.")

