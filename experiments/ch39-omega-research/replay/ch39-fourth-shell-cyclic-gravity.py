"""以三循环对称的第三步精确匹配 N+4 荷六的 K、W 和 D。"""
from pathlib import Path
from runpy import run_path

from flint import acb, acb_mat, arb, arb_mat


def verify(root):
    def lib(name):
        return run_path(str(root / ("ch39-" + name + ".py")))

    _, _, incoming, readout, (low, high, highest), weights, thermal = lib(
        "fourth-shell-two-gravity-rules")["shell_data"](root)
    moments, heat, field = thermal
    shell = lib("three-ray-chiral-boundary")
    sites = shell["geo"]["sites"]
    rotate = [sites.index((j,k,i)) for i,j,k in sites]
    rotation = acb_mat([[int(i == rotate[j]) for j in range(19)] for i in range(19)])
    # 目的谱简单且磁拉普拉斯与三循环对易，每条谱线的旋转本征值是一个三次单位根。
    adj = lambda a: a.conjugate().transpose()
    rotated = readout * rotation * adj(readout)
    omega = acb(-arb(1)/2, arb(3).sqrt()/2)
    groups = [[i for i in range(19) if rotated[i,i].contains(omega**(-k))] for k in range(3)]
    assert [len(g) for g in groups] == [7,6,6]
    assert sorted(i for group in groups for i in group) == list(range(19))
    assert all(rotated[i,j].contains(0) for i in range(19) for j in range(19) if i!=j)
    assert 15 in groups[0] and 17 in groups[0] and 18 in groups[1]

    overlap = readout * incoming * acb_mat([[(high-low).sqrt(),0],[0,(highest-low).sqrt()]])
    totals, maxima = [], []
    for group in groups:
        b = acb_mat([[overlap[i,j] for j in range(2)] for i in group])
        gram = adj(b)*b
        tr = gram.trace().real
        gap = ((gram[0,0].real-gram[1,1].real)**2 + 4*abs(gram[0,1])**2).sqrt()
        assert tr>gap>0
        totals.append(len(group)*low+tr)
        maxima.append(low+(tr+gap)/2)
    print("三个字符块的源迹与最大本征值", [v.str(14) for v in totals],
          [v.str(14) for v in maxima], flush=True)

    special = (15,17,18)
    bulk_weights = [sum(weights[i] for i in group if i not in special) for group in groups]
    w16,w18,w19 = [weights[i] for i in special]
    mz, mp = arb("91.1876"), arb("1.22089e19")
    p0,p1 = (mz/mp)**2,(arb("1.176397e14")/mp)**2
    weak = lib("native-rh-ray-pvm")["weak_series"](56)
    high_weak = lib("full-neutral-free-cutoff")["weak_column"]
    reference = arb("-.074256577194663365545")
    left,right = arb("1.4"), arb("2.1")
    def approximate(t):
        return w16*weak(t)+(1-w16)*reference
    assert approximate(left)>0>approximate(right)
    for _ in range(32):
        mid = (left+right)/2
        if approximate(mid)>0:
            left = mid
        else:
            assert approximate(mid)<0
            right = mid
    tbox = left.union(right)
    light = (p0*tbox).sqrt()
    assert light>low
    def fill(pair):
        return [(totals[0]-light-pair[0])/5, (totals[1]-pair[1])/5, totals[2]/6]

    def slope(a):
        ex = [(-z/a**2).exp() for _,_,z in field]
        return [a*sum(v*e for (_,v,_),e in zip(field,ex,strict=True))/arb.pi(),
                -a**3*sum(g*e for (g,_,_),e in zip(field,ex,strict=True))/(2*arb.pi())]

    def system(pair):
        bulk = fill(pair)
        rows = list(zip([*bulk_weights,w16,w18,w19],[*bulk,light,*pair],strict=True))
        values = [arb(-1),-arb("1e-122")]
        for w,a in rows:
            cols = moments(field,heat,a*a)
            for i in range(2):
                values[i] += w*cols[i]
        deriv = []
        for j,w in enumerate((w18,w19)):
            deriv.append([w*a-bulk_weights[j]*b/5
                          for a,b in zip(slope(pair[j]),slope(bulk[j]),strict=True)])
        return arb_mat([[z] for z in values]), arb_mat(deriv).transpose()

    center = [arb(".28"),arb(".60")]
    for _ in range(6):
        f,j = system(center)
        d = j.mid().solve(f.mid())
        center = [(a-b).mid() for a,b in zip(center,d.entries(),strict=True)]
    radius = arb("1e-10")
    pair = [z+arb(0,radius.upper()) for z in center]
    f,j = system(center)
    pre = j.mid().inv()
    assert not pre.det().contains(0)
    _,wide = system(pair)
    defect = arb_mat([[1,0],[0,1]])-pre*wide
    ratio = max((abs(z)/radius+sum(abs(v) for v in row)).upper()
                for z,row in zip((pre*f).entries(),defect.tolist(),strict=True))
    assert ratio<1,ratio

    bulk = fill(pair)
    for candidates,largest in zip(((light,pair[0],bulk[0]),(pair[1],bulk[1]),(bulk[2],)),
                                  maxima,strict=True):
        # 各块为 low I 加秩二正项。迹相同后，最大值条件与逐项下界蕴含全部主序条件。
        assert all(low<x<largest for x in candidates)
    fixed = sum(w*high_weak(z*z,p0,p1) for w,z in zip([*bulk_weights,w18,w19],[*bulk,*pair],strict=True))
    ends = [w16*weak(t)+fixed for t in (left,right)]
    assert ends[0]>0>ends[1]
    print("第16层能量 GeV、第18和19层振幅、三个字符的其余振幅",
          *[v.str(15) for v in (mz*tbox.sqrt(),*pair,*bulk)],flush=True)
    print("压缩比和弱差两端",ratio.str(8),*[v.str(10) for v in ends],flush=True)
    return pair,bulk,tbox


if __name__ == "__main__":
    verify(Path(__file__).resolve().parent)
