"""右外壳夸克中心要求输运将某个三循环不变态以超过一半的概率送出该子空间。"""
from pathlib import Path
from runpy import run_path

from flint import acb_mat, arb, ctx


def verify(root):
    model = run_path(str(root / "ch39-all-spaces-common-unitary-quark-right.py"))
    train = run_path(str(root / "ch39-common-embedding-central-right-ckm.py"))
    transport = run_path(str(root / "ch39-adjacent-pvm-transport.py"))
    shell = run_path(str(root / "ch39-three-ray-chiral-boundary.py"))
    source = 2_306_888_100_086
    _, labels, _ = train["verify_nodes"](root, ((source, ".0294"),), ".1", (),
                                          ((6,7,10,13), (6,7,10,12)))
    ctx.prec = 384
    adj = train["adj"]
    orbits = model["ORBITS"]
    basis = acb_mat([[1/arb(len(g)).sqrt() if i in g else 0 for g in orbits]
                     for i in range(19)])
    traces, characters = [], []
    for charge in (-4, 2):
        _, old = transport["spectrum"](shell, charge, source)
        _, new = transport["spectrum"](shell, charge, source+1)
        before = adj(basis)*old
        after = adj(basis)*new
        # 两端算子与三循环交换且谱单纯，各谱线在不变空间上的投影范数平方只能是 0 或 1。
        in_trivial = []
        other = []
        readout = []
        for j, weight in enumerate(labels[source,charge]):
            a = sum(abs(before[i,j])**2 for i in range(7))
            b = sum(abs(after[i,j])**2 for i in range(7))
            assert (a < arb(".01") or a > arb(".99")) and (b < arb(".01") or b > arb(".99"))
            (in_trivial if a > arb(".99") else other).append(weight.lower())
            if b > arb(".99"):
                readout.append((abs(after[2,j])**2 + abs(after[6,j])**2).lower())
        assert (len(in_trivial),len(other),len(readout)) == (7,12,7)
        traces.append(sum(a*b for a,b in zip(sorted(in_trivial), sorted(readout,reverse=True), strict=True)))
        characters.append(min(other))
    print("平凡块迹系数、非平凡字符系数",*[x.str(16) for x in traces+characters],flush=True)
    conditions = [(traces,arb(model["TRACE"])), (characters,arb(model["CHARACTER"]))]
    for (up,down),required in conditions:
        assert (1-arb(".50134"))*(up+25*down) > required
        p = arb(".23244")
        assert (1-2*p)**2*up + (1-p)*25*down > required
    print("两荷中较大的跨空间概率严格大于 0.50134；同一生成元的跨空间概率严格大于 0.23244。")


if __name__ == "__main__":
    verify(Path(__file__).resolve().parent)
