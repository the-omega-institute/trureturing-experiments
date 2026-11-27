"""固定四相位荷六立方根时，检验全部上型中心三平面的奇、底夸克质量。"""
from math import prod
from pathlib import Path
from runpy import run_path

from flint import acb, acb_mat, arb, ctx


def adj(a):
    return a.conjugate().transpose()


def line(v):
    return v*adj(v)*(1/(adj(v)*v)[0,0].real)


def frobenius(a):
    return sum(abs(z).upper()**2 for z in a.entries()).sqrt().upper()


def verify(root):
    witness = run_path(str(root/"ch39-rephased-involution-up-centers.py"))["verify"](root,return_matrices=True)
    base = run_path(str(root/"ch39-rank-two-generator-up-quark.py"))
    train = run_path(str(root/"ch39-common-embedding-central-right-ckm.py"))["verify_nodes"]
    shell = run_path(str(root/"ch39-three-ray-chiral-boundary.py"))
    spec = run_path(str(root/"ch39-adjacent-pvm-transport.py"))["spectrum"]
    n = base["N"] + 1
    _, labels, _ = train(root, ((n,".0294"),), ".1", (), ((6,7,10,13),(6,7,10,12)))
    ctx.prec = 384
    us, ut, ds, dt = [spec(shell,q,j)[1] for q,j in ((-4,n),(-4,n+1),(2,n),(2,n+1))]
    eye = base["eye"](19)
    # 用点值参考对合避免把相关参数在十九维矩阵乘积中当作独立误差累积。
    f = witness["first_center"]
    g = witness["source_differences"][1]
    gc = g.mid()
    g_gap = frobenius(g-gc)
    g_sine = g_gap/((adj(gc)*gc)[0,0].real.sqrt()-g_gap)
    f_sine = witness["first_motion"] / ((adj(f)*f)[0,0].real.sqrt()-witness["first_motion"])
    h0 = eye-2*(line(f)+line(gc))
    error = 2*(f_sine+g_sine)+witness["up_error"]
    def four(basis, weights):
        c = acb_mat([[basis[i,j] for j in range(19)] for i in (14,8,7,6)])
        return acb_mat([[c[i,j]*weights[j] for j in range(19)] for i in range(4)])*adj(c)
    uq = witness["up_center"]
    down = four(dt,base["read"](ds,dt,labels[n,2],h0*uq))
    up = four(ut,base["read"](us,ut,labels[n,-4],uq))
    roots = [arb(1), acb(-1,arb(3).sqrt())/2, acb(-1,-arb(3).sqrt())/2]
    axes = [acb_mat([[0],[1],[z],[z*z]])*(1/arb(3).sqrt()) for z in roots[1:]]
    down_factor = max(w.upper() for w in labels[n,2]).sqrt()
    up_factor = max(w.upper() for w in labels[n,-4]).sqrt()
    def bound(v, radius):
        return (max(arb(0),v.sqrt()-radius)**2).union((v.sqrt()+radius)**2)
    dvals = sorted((z.real for z in down.eig(algorithm="rump")),key=lambda z:z.mid())
    down_boxes = [bound(v,down_factor*error) for v in dvals]
    levels = witness["levels"]
    masses = list(map(arb,base["CENTERS"]))
    for axis in axes:
        uc = bound((adj(axis)*up*axis)[0,0].real,up_factor*witness["up_error"])
        dc = bound((adj(axis)*down*axis)[0,0].real,down_factor*error)
        if dc > arb(".1"):
            assert uc.overlaps(levels[1]) and all(not uc.overlaps(levels[j]) for j in (0,2,3))
            a = levels[1]
            prob = prod(a-m for m in masses)/prod(a-v for j,v in enumerate(levels) if j!=1)
            bottom = dc*(1-prob)
            break
    else:
        raise AssertionError("未能辨认高下型字符")
    assert down_boxes[1] > arb(".0045") and bottom > arb(".118")
    print("固定根的下型四站谱外包",[(x.lower().str(9),x.upper().str(9)) for x in down_boxes])
    print("高下型字符在上型中心法向上的概率及底质量下界",prob.str(8),bottom.lower().str(12))
    return down_boxes, bottom


if __name__ == "__main__":
    verify(Path(__file__).resolve().parent)
