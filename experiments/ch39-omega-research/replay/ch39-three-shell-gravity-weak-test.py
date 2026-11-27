"""按三壳引力六式的同一角盒，检验三个校准壳在匹配点的弱差。"""
from pathlib import Path
from runpy import run_path

from flint import arb


def experiment(root, fifth, columns, levels):
    def lib(name):
        return run_path(str(root / ("ch39-" + name + ".py")))

    tri = lib("three-shell-gravity-radius56")
    box = tri["experiment"](root, fifth, columns, levels)
    ordered = lib("two-shell-ordered-bond-gates")
    transport, *_ = ordered["transport"](fifth*columns, box, parameters=tri["PARAMETERS"],
                                        active=tri["ACTIVE"], layers=ordered["LAYERS"] + tri["EXTRA_LAYERS"])
    low, high, top = levels
    weak = lib("full-neutral-free-cutoff")["weak_column"]
    p0 = (arb("91.1876")/arb("1.22089e19"))**2
    p1 = (arb("1.176397e14")/arb("1.22089e19"))**2
    reference = arb("-.074256577194663365546").union(arb("-.074256577194663365544"))
    for number, displacement in zip(range(2_306_888_100_092, 2_306_888_100_095),
                                    (".00000822", ".00001266", ".00001839"), strict=True):
        _, _, _, basis, _, weights, _, _ = lib("fifth-seventh-shortest-site-gates")["geometry"](
            root, fifth, columns, levels, number=number
        )
        projection = basis.conjugate().transpose()*transport
        cutoff = []
        for i in range(19):
            modes = [abs(projection[i,j]) for j in range(2)]
            left = low.lower()+sum(g.lower()*max(arb(0), mode.lower())**2
                                   for g,mode in zip((high-low,top-low), modes, strict=True))
            cutoff.append(left**2)
        print("最小平方截止与谱权和", min(cutoff).str(12), sum(weights).str(8), flush=True)
        assert min(cutoff) > p0 and sum(weights).contains(1)
        # 谱权的和精确为一；按共同参考值居中后，余项按谱权绝对值包围。
        discrepancy = reference + sum(w*(weak(x,p0,p1)-reference) for w,x in zip(weights,cutoff,strict=True))
        print("壳、弱差的有向区间和最小平方截止下界", number, discrepancy.lower().str(18),
              discrepancy.upper().str(18), min(cutoff).str(12), flush=True)
        assert arb("-.07425657719466338") < discrepancy < arb("-.07425657719466335")
        weight_norm = sum(w.abs_upper() for w in weights)
        necessary = (arb("7.468e-36")/(reference.abs_lower()/weight_norm-arb("2.398e-21"))).sqrt()
        cap = max(p0.sqrt().upper(), necessary.upper())
        distance = min(cutoff).sqrt().lower()-cap
        assert distance > arb(displacement)
        print("弱差归零所需的一条谱线的最高振幅及与现解的算子范数距离下界",
              cap.str(14), distance.str(14), flush=True)


if __name__ == "__main__":
    root = Path(__file__).resolve().parent
    source = run_path(str(root / "ch39-fifth-full-positive-joint.py"))
    experiment(root, *source["verify"](root, return_operator="spectral"))
