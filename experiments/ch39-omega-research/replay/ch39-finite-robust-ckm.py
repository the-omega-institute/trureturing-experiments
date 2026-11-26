"""Finite CKM and quark hierarchy on all diagonal choices, in orthonormal ring labels."""
from pathlib import Path
import runpy
from flint import arb

finite=runpy.run_path(str(Path(__file__).with_name("ch39-finite-orthogonal-ckm.py")))
fit,common=finite["fit"],finite["common"]
weights=(
    (".005542097281",".015787288971",".164285816006",".298811468278",
     ".963799203568",".340719042261",".340734221108","1.08672805250",".011775145159",".537620773878"),
    (".080003313538",".001722155453",".002628511184",".119258420177",
     "2.22830530721",".096599394566",".877919639392",".004413973180",".349545824757",".362256612626"))


def matrices(model,weight):
    w=list(map(arb,weight))
    k=fit["evaluate"](model,w[:4],(0,0))[2]
    for end in (1,2):
        for j in range(3):
            finite["add_edge"](k,fit["inner"][j],fit["outer"][(j+end)%3],w[4+3*(end-1)+j])
    inn=finite["fix_phase"](k,(7,6,8),arb.pi()-arb(".002"))
    base=finite["fix_phase"](k,(11,12,18),arb(".03606860"))
    degrees,margin=finite["flux"]["degrees"](k)
    assert degrees==(3,[3]) and margin>arb(".0001")
    t=-fit["block"](k,fit["ii"],fit["oo"])
    print("inner/bottom strengths and angular cut margin",inn.str(10),base.str(10),margin.str(8))
    return t*t


if __name__=="__main__":
    model=fit["model"](*common["angles"])
    vals=finite["observables"](*(matrices(model,w) for w in weights))
    goals=common["up"]+list(reversed(common["down"]))+list(map(arb,(".22501",".04183",".003732","1.147")))
    assert all(abs(v/target-1)<arb("1e-6") for v,target in zip(vals,goals))
    assert arb(".0000300")<vals[-1]<arb(".0000325")
    print("mass ratios, mixing parameters and J",*(v.str(11) for v in vals))
