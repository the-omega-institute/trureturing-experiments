"""冻结中央零点的正振幅谱权，在两个相邻实际编号重算六个 Yukawa 与混合。"""
from pathlib import Path
from runpy import run_path

from flint import acb_mat, arb, ctx


def verify(root):
    native = run_path(str(root / "ch39-free-four-site-native-positive-spectrum.py"))
    ckm = run_path(str(root / "ch39-shell-clock-ckm-blocks.py"))
    spectrum = run_path(str(root / "ch39-three-gram-source-flavor.py"))["projections"]
    ctx.prec = 320
    center = 2_306_888_100_085
    rows_at_center = native["native"](root, center)
    adj = native["adj"]
    v = acb_mat(ckm["standard_matrix"](*[arb(x).sin() for x in (".22704",".04275",".00377")], arb("1.139")))
    mass = [list(map(arb, seq)) for seq in ((".00000615",".00311",".8616"),(".0000135",".000268",".01401"))]
    width = [list(map(arb, seq)) for seq in ((".00000014",".00005",".0043"),(".0000002",".000003",".00009"))]
    a, angle = arb(".02145737"), arb("-.44525")
    frame = native["physical_frame"](a, angle)
    up, down = [acb_mat([[y[i]/y[2] if i==j else 0 for j in range(3)] for i in range(3)]) for y in mass]
    up = adj(v)*up*v
    common = native["target_basis"](up, down)
    responses = []
    for rows, normalized, indices in zip(rows_at_center, (up,down), ((1,6,13,16),(1,6,10,16)), strict=True):
        target = adj(common)*normalized*common
        _, weights = native["decompose"](rows, target, a, angle, indices, return_weights=True)
        responses.append(weights)

    for offset in (-1,1):
        projectors = []
        missed = 0
        for rows, weights, targets, errors in zip(native["native"](root,center+offset), responses, mass, width, strict=True):
            block = sum((w*r for w,r in zip(weights,rows,strict=True)), acb_mat(4,4))
            amplitude = adj(frame)*block*frame*targets[2]
            levels, modes = spectrum(amplitude*amplitude)
            singular = [x.sqrt() for x in levels]
            flags = [abs(x-y)>err for x,y,err in zip(singular,targets,errors,strict=True)]
            assert all(flag or abs(x-y)<err for flag,x,y,err in zip(flags,singular,targets,errors,strict=True))
            missed += sum(flags)
            print("偏移、三个汤川耦合、逐项错窗",offset,*[x.str(9) for x in singular], flags)
            projectors.append(modes)
        pu,pd = projectors
        vub = (pu[0]*pd[2]).trace().real.sqrt()
        j = (pu[0]*pd[0]*pu[1]*pd[1]).trace().imag
        print("错窗质量数、Vub、带符号 J",missed,vub.str(10),j.str(10))
        s12 = ((pu[0]*pd[1]).trace().real / (1-vub*vub)).sqrt()
        s23 = ((pu[1]*pd[2]).trace().real / (1-vub*vub)).sqrt()
        print("三个角",s12.asin().str(11),s23.asin().str(11),vub.asin().str(11))
        assert missed == (3 if offset == -1 else 4)
        assert abs(s12.asin()-arb(".22704")) > arb(".00082")


if __name__ == "__main__":
    verify(Path(__file__).resolve().parent)
