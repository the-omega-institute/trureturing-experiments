"""混合顶点相消需要偏离三类原生面数的份额。整数数据单位是 10^-10。"""
import json
from pathlib import Path
from runpy import run_path

from flint import acb, arb, ctx


def verify():
    root = Path(__file__).resolve().parent
    def load(name):
        return run_path(str(root / ("ch39-"+name+".py")))

    fixed = load("height-independent-flavor-vacuum")
    lib = load("mincut-global-winding")
    shell = load("three-ray-chiral-boundary")
    flavor = load("chiral-neutral-flavor-phase")
    angles = load("mincut-reference-ckm-action")
    ctx.prec = 224
    tri = lib["terms"](shell)
    x = arb("-.01046")
    y = x-angles["reference_angles"](flavor)[0]
    h = acb.zeta_zero(2_306_888_100_085).imag
    phase = shell["geo"]["phases"](h)
    actual = (phase[724].arg(), phase[804].arg())
    # 对偶不等式使用单位辐边下的全局能序；顶点则取实际的软边亏额。
    delta = arb(23)/25*arb(2)**-20/(h*h+shell["soft"]["point"]["study"]["offset"])
    bare = {q:lib["coeff"](tri,q,arb(0)) for q in range(1,106)}
    spoke = {q:lib["coeff"](tri,q,delta) for q in bare}
    gradients = {q:fixed["gradient_columns"](lib,tri,q,x,y) for q in bare}
    data = json.loads((root/"ch39-proton-class-share-cuts.json").read_text())
    scale = 10**10
    vertices = {(scale,0,0), (0,scale,0), (2550000000,0,7450000000),
                (0,2550000000,7450000000)}
    assert {tuple(row["class"]) for row in data} == vertices
    for row in data:
        shares = [arb(i)/scale for i in row["class"]]
        equality = [arb(i)/scale for i in row["eq"]]
        cuts = [[arb(i)/scale for i in point] for point in row["cuts"]]
        assert all(w>0 for _,_,w in cuts)
        for a,b in ((x,y),actual):
            lower = []
            for q,cs in bare.items():
                f = sum(s*(z*acb(0,q*t).exp()).real/n for s,z,t,n in
                        zip(shares,spoke[q],(a,b,a-b),(42,42,12),strict=True))
                at_vacuum = lib["energy"](q,cs,x,y)
                compares = sum(w*(lib["energy"](q,cs,i,j)-at_vacuum) for i,j,w in cuts)
                correction = sum(w*g for w,g in zip(equality,gradients[q],strict=True))
                lower.append(f-compares-correction)
            bound = min(v.lower() for v in lower)
            print("整数顶点、逐阶下界",row["class"],bound.str(14))
            assert bound > arb(".00009")


if __name__ == "__main__":
    verify()
