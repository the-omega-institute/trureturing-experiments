"""在同一实际背景上检验每条原边的整个非负电导半轴。"""
from pathlib import Path
from runpy import run_path
import sys

from flint import acb, acb_mat, arb, ctx


def rank_one_value(inverse, phase, i, j):
    # b=e_i-conjugate(phase)e_j，故 b*Gb=Gii+Gjj-2Re(conjugate(phase)Gij)。
    return inverse[i,i].real + inverse[j,j].real - 2*(phase.conjugate()*inverse[i,j]).real


def certify(native, edges, levels, eye, radius):
    lower_cut, upper_cut = levels[2]+radius, levels[3]-radius
    assert lower_cut < levels[3] and upper_cut > levels[2]
    lower = (native-lower_cut*eye).inv()
    upper = (native-upper_cut*eye).inv()
    small, large = [], []
    for i,j,_,_ in edges:
        weight = abs(native[i,j])
        assert weight > 0
        phase = -native[i,j]/weight
        small.append(1 - weight*rank_one_value(lower,phase,i,j))
        large.append(rank_one_value(upper,phase,i,j))
    return min(z.lower() for z in small), min(z.lower() for z in large)


def verify(root):
    shell = run_path(str(root/"ch39-three-ray-chiral-boundary.py"))
    ctx.prec = 192
    h = acb.zeta_zero(22).imag
    edges = shell["geo"]["edges"]
    assert len(edges) == 135
    eye = acb_mat([[int(i==j) for j in range(19)] for i in range(19)])
    for q in (1,-4,2,-3,6):
        native = shell["laplacian"](q,h)
        raw = native.eig(multiple=True)
        levels = sorted([z.real for z in raw],key=float)
        assert len(levels)==19 and all(z.imag.contains(0) for z in raw)
        assert all(v<w for v,w in zip(levels,levels[1:]))
        for num in range(10,0,-1):
            radius = arb(num)/10
            if radius >= levels[3]-levels[2]:
                continue
            decay, grow = certify(native,edges,levels,eye,radius)
            if decay>0 and grow>0:
                break
        else:
            raise AssertionError("未在本荷的所有边取得统一正余量")
        print("第 22 枚，荷、135 条单边全权谱隙严格下界、减权和增权两判据最小余量",
              q,radius.str(8),decay.str(12),grow.str(12))


if __name__=="__main__":
    root = Path(sys.argv[1]).resolve() if len(sys.argv)>1 else Path(__file__).resolve().parent
    del sys.argv[1:]
    verify(root)
