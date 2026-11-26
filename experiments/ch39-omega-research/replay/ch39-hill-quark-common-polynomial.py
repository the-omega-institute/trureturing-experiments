"""三十次实多项式的十九个荷六截止条件及上下夸克的十二个压缩条件。"""
from pathlib import Path
from runpy import run_path

from flint import acb, acb_mat, arb, arb_mat, ctx


def product(factors):
    value = arb(1)
    for term in factors:
        value *= term
    return value


def four_coordinates(block):
    spoke, edge = block[0,1], block[1,2]
    return [block[0,0].real, block[1,1].real, spoke.real, spoke.imag, edge.real, edge.imag]


def compression(frame, common):
    columns = []
    for i in range(6):
        p, r, a, b, c, d = [int(i==j) for j in range(6)]
        z, w = acb(a,b), acb(c,d)
        block = acb_mat([[p,z,z,z],[z.conjugate(),r,w,w.conjugate()],
                         [z.conjugate(),w.conjugate(),r,w],[z.conjugate(),w,w.conjugate(),r]])
        columns.append(common["coordinates"](common["adj"](frame)*block*frame))
    return arb_mat(columns).transpose()


def verify(root):
    ray = run_path(str(root / "ch39-three-ray-chiral-boundary.py"))
    common = run_path(str(root / "ch39-signed-three-source-quark-pair.py"))
    native = run_path(str(root / "ch39-native-rh-ray-pvm.py"))
    ckm = run_path(str(root / "ch39-shell-clock-ckm-blocks.py"))
    flavor = run_path(str(root / "ch39-three-gram-source-flavor.py"))
    leakage = run_path(str(root / "ch39-three-source-amplitude-leakage.py"))
    run_path(str(root / "ch39-charge-six-hill-cutoff.py"))["verify"](root)
    ctx.prec = 416
    height = acb.zeta_zero(2_306_888_100_085).imag
    pos = [ray["soft"]["pole"], *ray["boundary_sites"]()]
    assert pos == [14, 8, 7, 6]
    levels = ray["laplacian"](6, height).eig(algorithm="rump")
    assert all(z.imag.contains(0) for z in levels)
    nodes = sorted(((v.real - 18)/8 for v in levels), key=lambda v: v.mid())
    assert all(a < b for a, b in zip(nodes, nodes[1:]))
    slopes = [product(x-y for j,y in enumerate(nodes) if j!=i) for i,x in enumerate(nodes)]
    frame = common["adapted_frame"](native["gram_cholesky"](root))
    adj = common["adj"]
    rank_rows = []
    lagrange_rows, spectral_rows = [], []
    for q in (-4, 2):
        matrix = ray["laplacian"](q, height)
        vals, vectors = matrix.eig(right=True, algorithm="rump")
        assert all(v.imag.contains(0) for v in vals)
        weighted = []
        for i, level in enumerate(vals):
            t = (level.real-18)/8
            chi = product(t-a for a in nodes)
            assert not chi.contains(0)
            norm = sum(abs(vectors[k,i])**2 for k in range(19))
            v = acb_mat([[vectors[k,i]] for k in pos])
            proj = (v*adj(v))*(1/norm)
            cardinal = [chi/((t-x)*slope) for x,slope in zip(nodes,slopes)]
            weighted.append((t, chi, proj, cardinal))
        rows = []
        for degree in range(12):
            block = sum((p*(chi*t**degree) for t,chi,p,_ in weighted), acb_mat(4,4))
            rows.append(four_coordinates(block))
        rank_rows.extend(arb_mat(rows).transpose().tolist())
        columns = []
        for j in range(19):
            block = sum((p*cardinal[j] for _,_,p,cardinal in weighted), acb_mat(4,4))
            columns.append(four_coordinates(block))
        lagrange_rows.extend(arb_mat(columns).transpose().tolist())
        spectral_rows.append(weighted)
    matrix = arb_mat(rank_rows)
    mid = arb_mat([[v.mid() for v in row] for row in rank_rows])
    identity = arb_mat([[int(i==j) for j in range(12)] for i in range(12)])
    defect = mid.inv()*matrix - identity
    norm = max(sum(abs(x).upper() for x in row) for row in defect.tolist())
    print("四站十二阶矩阵的中点行列式和预条件缺陷",mid.det().str(14),norm.str(14))
    assert norm < arb("4.16e-20")
    change = compression(frame, common)
    assert not change.det().contains(0)
    print("从四站六坐标到适配三源标架的行列式",change.det().str(12))
    h = arb("5.66600e-35").union(arb("5.66604e-35"))
    kappa = arb(".635595244").union(arb(".635595246"))
    s = arb("3.291030650e-40").union(arb("3.291030654e-40"))
    powers = [((8*t+18)/(8*nodes[-1]+18))**6000 for t in nodes[:-1]] + [arb(1)]
    sites = [h+kappa*u/(u+s) for u in powers]
    assert all(0<x<1 for x in sites)
    targets = common["target_blocks"](ckm,flavor)
    amplitudes, up, down, v = leakage["target_amplitudes"](ckm,common)
    left = []
    for a,masses in zip(amplitudes,(up,down),strict=True):
        sq, projectors = flavor["projections"](a*a)
        assert all(l>0 and l.overlaps(y**2) for l,y in zip(sq,masses,strict=True))
        left.append(projectors)
    pu,pd = left
    assert all((pu[i]*pd[j]).trace().real.overlaps(abs(v[i,j])**2) for i in range(3) for j in range(3))
    assert (pu[0]*pd[0]*pu[1]*pd[1]).trace().imag.overlaps(
        (v[0,0]*v[1,1]*v[0,1].conjugate()*v[1,0].conjugate()).imag)
    lagrange = arb_mat(lagrange_rows)
    for name, cutoff, blocks in (("平方质量",sites,targets), ("振幅",[x.sqrt() for x in sites],amplitudes)):
        rhs = []
        for block in blocks:
            target = arb_mat([[x] for x in common["coordinates"](block)])
            rhs.extend(change.solve(target,algorithm="precond").tolist())
        b = arb_mat(rhs)
        cut = arb_mat([[v] for v in cutoff])
        def functional(monomials, cardinal):
            # 先对目标做线性消元，再代入 Hill 的区间，保留同一输入在各层间的相消。
            dual = matrix.transpose().solve(monomials,algorithm="precond").transpose()
            return (dual*b + (cardinal-dual*lagrange)*cut)[0,0]
        correction = [functional(arb_mat([[int(i==j)] for i in range(12)]),arb_mat(1,19)) for j in range(12)]
        leading = correction[-1]
        if name=="振幅":
            assert arb("-1.13089e14") < leading < arb("-1.13088e14")
        else:
            assert arb("6.89788e13") < leading < arb("6.89790e13")
        print(name,"χ6 乘子的十二个升幂系数", *(v.str(10) for v in correction))
        for q,data,block in zip((-4,2),spectral_rows,blocks,strict=True):
            weights = [functional(arb_mat([[chi*t**j] for j in range(12)]),arb_mat([cardinal]))
                       for t,chi,p,cardinal in data]
            counts = sum(x>0 for x in weights),sum(x<0 for x in weights)
            expected = (7,12) if name=="振幅" else ((11,8) if q==-4 else (8,11))
            assert counts == expected
            print(name,q,"严格为正、负的谱点数",*counts)
            if name=="振幅":
                square = sum((p*(w*w) for (_,_,p,_),w in zip(data,weights)),acb_mat(4,4))
                leaked = (adj(frame)*square*frame-block*block).trace().real
                bounds = ("223028.94","223028.96") if q==-4 else ("1819.3037","1819.3039")
                assert arb(bounds[0])<leaked<arb(bounds[1])
                print("先平方再压缩多出的迹",q,leaked.str(10))


if __name__ == "__main__":
    verify(Path(__file__).resolve().parent)
