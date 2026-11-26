"""固定正目标味基时，十九层实汤川振幅的最小最大模。"""
from pathlib import Path
from runpy import run_path

from flint import acb, arb, arb_mat, ctx


def groups(charge):
    plus = (2,4,8,11,14,17)
    minus = (3,5,9,12,15,19) if charge == -4 else (3,5,9,13,15,19)
    return [tuple(j-1 for j in g) for g in (plus, minus)] + [
        (j-1,) for j in range(1,20) if j not in plus+minus]


def certify(target, data, groups_in_order):
    # 每个非平凡三循环特征空间在四站上只有一个方向。用 yz 虚部的符号区分。
    ray_signs = []
    for indexes in groups_in_order[:2]:
        yz = [data[j][1][5,0] for j in indexes]
        assert all(data[j][1][1,0] > 0 for j in indexes)
        sign = 1 if yz[0] > 0 else -1
        assert all(sign*v > 0 for v in yz)
        ray_signs.append(sign)
    assert ray_signs[0] == -ray_signs[1]
    columns = [sum((data[i][1] for i in group), arb_mat(6,1)) for group in groups_in_order]
    free = (0,1,2,7,8)
    active = (3,4,5,6)
    saturation = -sum((columns[i] for i in active), arb_mat(6,1))
    system = arb_mat([columns[i].entries() for i in free] + [saturation.entries()]).transpose()
    assert not system.det().contains(0)
    result = system.solve(target, algorithm="precond").entries()
    peak = result[-1]
    assert all(abs(x) < peak for x in result[:-1])
    all_weights = [-peak] * 19
    for group, value in zip(free,result[:-1],strict=True):
        for i in groups_in_order[group]:
            all_weights[i] = value
    residual = sum((w*d for w,(_,d) in zip(all_weights,data,strict=True)), arb_mat(6,1))-target
    assert residual.contains(arb_mat(6,1))

    e6 = arb_mat([[int(i == 5)] for i in range(6)])
    dual = system.transpose().solve(e6, algorithm="precond")
    support = [(dual.transpose()*d)[0,0] for _,d in data]
    assert all(support[j] < 0 for i in active for j in groups_in_order[i])
    lower = (dual.transpose()*target)[0,0]/sum(abs(v) for v in support)
    assert lower.overlaps(peak)
    print("自由组值、最小峰值",*(v.str(16) for v in result))
    print("对偶向量",*(v.str(12) for v in dual.entries()))
    print("全十九变量对偶下界",lower.str(16))
    return lower, peak, all_weights


def verify(root):
    ray=run_path(str(root/"ch39-three-ray-chiral-boundary.py"))
    gram=run_path(str(root/"ch39-native-rh-ray-pvm.py"))
    ckm=run_path(str(root/"ch39-shell-clock-ckm-blocks.py"))
    common=run_path(str(root/"ch39-signed-three-source-quark-pair.py"))
    leak=run_path(str(root/"ch39-three-source-amplitude-leakage.py"))
    ctx.prec=320
    height=acb.zeta_zero(2_306_888_100_085).imag
    frame=common["adapted_frame"](gram["gram_cholesky"](root))
    sites=[ray["soft"]["pole"],*ray["boundary_sites"]()]
    targets,*_=leak["target_amplitudes"](ckm,common)
    six = ray["laplacian"](6,height).eig(algorithm="rump")
    assert all(v.imag.contains(0) for v in six)
    six = sorted((v.real for v in six),key=lambda v:v.mid())
    h = arb("5.66600e-35").union(arb("5.66604e-35"))
    kappa = arb(".635595244").union(arb(".635595246"))
    s = arb("3.291030650e-40").union(arb("3.291030654e-40"))
    powers = [(v/six[-1])**6000 for v in six[:-1]]+[arb(1)]
    cutoffs = [h+kappa*r/(r+s) for r in powers]
    assert all(0<x<1 for x in cutoffs)
    interpolation = [(point,x.sqrt()) for point,x in zip(six,cutoffs,strict=True)]
    for q,target in zip((-4,2),targets):
        levels,_,data,_=leak["moment_frame"](ray["laplacian"](q,height),frame,sites,common)
        print("荷、活跃层序号",q,[groups(q)[i][0]+1 for i in (3,4,5,6)])
        lower, peak, weights = certify(leak["real_coordinates"](target,common), data, groups(q))
        lo,hi = ("43.1079147041755","43.1079147041757") if q==-4 else ("2.64783908069605","2.64783908069609")
        assert arb(lo)<lower and peak<arb(hi)
        print("超出单位幅度上限", (peak-1).str(15))
        interpolation.extend(zip(levels,weights,strict=True))
    interpolation.sort(key=lambda pair:pair[0].mid())
    nodes = [point for point,value in interpolation]
    assert len(nodes)==57 and all(a < b for a,b in zip(nodes,nodes[1:]))
    print("57 个跨荷插值节点的最小间距",min((b-a).lower() for a,b in zip(nodes,nodes[1:])).str(12))
    return interpolation


if __name__=="__main__":
    verify(Path(__file__).resolve().parent)
