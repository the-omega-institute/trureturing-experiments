"""有界的同一正有理谱读出与两枚相邻零点的荷六截止。"""
from pathlib import Path
from runpy import run_path

from flint import acb, acb_mat, arb, ctx


def ordered(shell, charge, height):
    values = shell["laplacian"](charge, height).eig(algorithm="rump")
    assert all(v.imag.contains(0) for v in values)
    levels = sorted([v.real for v in values], key=lambda x: x.mid())
    assert all(x < y for x,y in zip(levels,levels[1:]))
    return levels


def verify(root):
    native = run_path(str(root / "ch39-free-four-site-native-positive-spectrum.py"))
    ckm = run_path(str(root / "ch39-shell-clock-ckm-blocks.py"))
    shell = run_path(str(root / "ch39-three-ray-chiral-boundary.py"))
    ctx.prec = 320
    center = 2_306_888_100_085
    height = acb.zeta_zero(center).imag
    adj = native["adj"]
    v = acb_mat(ckm["standard_matrix"](*[arb(x).sin() for x in (".22704",".04275",".00377")],arb("1.139")))
    sector_masses = [list(map(arb, row)) for row in ((".00000615",".00311",".8616"),(".0000135",".000268",".01401"))]
    up,down = [acb_mat([[y[i]/y[2] if i==j else 0 for j in range(3)] for i in range(3)]) for y in sector_masses]
    up = adj(v)*up*v
    q = native["target_basis"](up,down)
    interpolation = []
    for charges,rows,m,ys,indices in zip((-4,2), native["native"](root,center), (up,down), sector_masses,
                                          ((1,6,13,16),(1,6,10,16)),strict=True):
        _,weights = native["decompose"](rows,adj(q)*m*q,arb(".02145737"),arb("-.44525"),indices,return_weights=True)
        interpolation.extend(zip(ordered(shell,charges,height),[w*ys[-1] for w in weights],strict=True))
    six = ordered(shell,6,height)
    h = arb("5.66600e-35").union(arb("5.66604e-35"))
    cap = arb(".635595244").union(arb(".635595246"))
    s = arb("3.291030650e-40").union(arb("3.291030654e-40"))
    for i,t in enumerate(six):
        ratio = arb(1) if i == 18 else (t/six[-1])**6000
        value = (h+cap*ratio/(ratio+s)).sqrt()
        assert 0 < value < 1
        interpolation.append((t,value))
    interpolation.sort(key=lambda pair:pair[0].mid())
    assert len(interpolation)==57 and all(0<w<arb("1.77") for t,w in interpolation)
    gaps = [v[0]-u[0] for u,v in zip(interpolation,interpolation[1:])]
    assert min(g.lower() for g in gaps) > arb(".000205")
    # 带宽取固定二进分数。
    width = arb(2)**-48
    background = arb(".5")
    for shift in (-1,1):
        next_height = acb.zeta_zero(center+shift).imag
        values, near = [], []
        for charge in (6,-4,2):
            for x in ordered(shell,charge,next_height):
                distances = [abs(x-t) for t,w in interpolation]
                near.append(min(d.lower() for d in distances))
                assert near[-1] > 0
                weights = [(width/d)**2 for d in distances]
                response = (background+sum(v*w for v,(_,w) in zip(weights,interpolation,strict=True)))/(1+sum(weights))
                values.append(response)
        spread = max(abs(v-background).upper() for v in values)
        print("偏移、最近中央节点间距、振幅偏离 1/2 的最大包络",shift,min(near).str(10),spread.str(10))
        distance_floor, error_cap = (("8.28e-7", "9.2e-18") if shift == -1 else ("3.35e-6", "5.7e-19"))
        assert min(near) > arb(distance_floor) and spread < arb(error_cap)
    print("中央 57 项的严格正下界及上界", min(v.lower() for t,v in interpolation).str(10),
          max(v.upper() for t,v in interpolation).str(10))


if __name__ == "__main__":
    verify(Path(__file__).resolve().parent)
