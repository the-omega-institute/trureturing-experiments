"""固定第五壳源，七壳弱差同时为零要求酉输运至少移动八个站位。"""
from itertools import combinations, product
from pathlib import Path
from runpy import run_path

from flint import acb_mat, arb, ctx


def adj(matrix):
    return matrix.conjugate().transpose()


def norm(row):
    return sum(abs(z).upper()**2 for z in row).sqrt().upper()


def restrict(matrix, indices):
    return acb_mat([matrix.tolist()[i] for i in indices])


def smaller_eigenvalue(matrix):
    a, b, d = matrix[0,0].real, matrix[0,1], matrix[1,1].real
    return (a+d-((a-d)**2+4*abs(b)**2).sqrt())/2


def experiment(root, fifth, columns, levels):
    def lib(name):
        return run_path(str(root/("ch39-"+name+".py")))

    graph = lib("three-ray-chiral-boundary")
    frame = lib("cyclic-two-characters-nine")["characters"](graph["geo"]["sites"])[0][0]
    orbits = [[i for i in range(19) if abs(frame[i,k]) > arb(".5")] for k in range(7)]
    choices = lib("seven-cyclic-weak-horizon")["CANDIDATES"]
    floor, heavy, top = levels
    ctx.prec = 256
    weighted = fifth*columns*acb_mat([[(heavy-floor).sqrt(),0],[0,(top-floor).sqrt()]])
    assert sorted(map(len,orbits)) == [1,3,3,3,3,3,3]
    assert sorted(i for orbit in orbits for i in orbit) == list(range(19))
    source = adj(frame)*weighted
    rows = []
    for number in range(2_306_888_100_092,2_306_888_100_099):
        _, basis = lib("adjacent-pvm-transport")["spectrum"](graph, 6, number)
        rows.append((adj(basis)*frame).tolist())
    first = [shell[0] for shell in rows]
    radius = (arb(3000)/arb("1.22089e19")-floor).sqrt()
    total = 0
    # 六条被触及轨道中，一条给两个站，其余各给一个。固定轨道决定中心。
    for absent, drop, bound, expected in [(3,3,".15",1215),(4,2,".15",1215),(6,2,".100",4374)]:
        active = [k for k in range(7) if k != absent]
        lines = [row for j,row in enumerate(first) if j != drop]
        inverse = acb_mat([[row[k] for k in active] for row in lines]).inv()
        fixed = acb_mat([[row[absent]] for row in lines])*restrict(source,[absent])
        center = -inverse*fixed
        radii = [radius*sum(abs(z) for z in row) for row in inverse.tolist()]
        count, gap = 0, arb(10)
        for double in (k for k in active if len(orbits[k]) > 1):
            groups = [list(combinations(orbits[k],2 if k==double else 1)) for k in active]
            scale = [(arb(len(orbits[k]))/len(group[0])).sqrt() for k,group in zip(active,groups,strict=True)]
            bounds = [a*b for a,b in zip(scale,radii,strict=True)]
            for support in product(*groups):
                local = restrict(weighted,[site for group in support for site in group])
                predicted = []
                for j,(k,group) in enumerate(zip(active,support,strict=True)):
                    predicted.append([sum(weighted[i,h] for i in group)/arb(len(group)).sqrt()
                                      +scale[j]*(center[j,h]-source[k,h]) for h in (0,1)])
                # 每个轨道上只保留被选站的归一化平均行。压缩得到的 Gram 须不大于原局部 Gram。
                y = acb_mat(predicted)
                delta = adj(local)*local-adj(y)*y
                error = sum((2*norm(row)+d)*d for row,d in zip(predicted,bounds,strict=True))
                value = (-smaller_eigenvalue(delta)-error).lower()
                assert value.is_finite() and value > arb(bound), (support,value)
                gap = min(gap,value)
                count += 1
        assert count == expected and gap > arb(bound), (count,gap)
        total += count
        print("不动轨道、七站支撑数、压缩 Gram 的最小违反量",orbits[absent],count,gap.str(12),flush=True)
    assert total == 6804

    # 七条轨道各一站，混合选线的每个逆矩阵均给同一中心 0 周围的行范数包。
    sizes = [arb(len(orbit)).sqrt() for orbit in orbits]
    caps = [arb(0) for _ in orbits]
    count = 0
    for ids in product(*choices):
        if any(ids):
            inverse = acb_mat([r[i] for r,i in zip(rows,ids,strict=True)]).inv()
            for j,row in enumerate(inverse.tolist()):
                value = (radius*sizes[j]*sum(abs(z) for z in row)).upper()
                assert value.is_finite()
                caps[j] = max(caps[j],value)
            count += 1
    assert count == 971
    # 全第一线的病态方向：释放第七个轨道均值，由六壳求其他均值，再正交消去这一自由方向。
    lines = [row for j,row in enumerate(first) if j!=2]
    inverse = acb_mat([row[:6] for row in lines]).inv()
    ratios = (-inverse*acb_mat([[row[6]] for row in lines])).entries()+[arb(1)]
    direction = acb_mat([[a*b] for a,b in zip(sizes,ratios,strict=True)])
    projector = direction*adj(direction)/(adj(direction)*direction)[0,0]
    first_error = norm([radius*s*sum(abs(z) for z in row)
                        for s,row in zip(sizes[:6],inverse.tolist(),strict=True)])
    mixed_gap, first_gap, count = arb(10), arb(10), 0
    for support in product(*orbits):
        local = restrict(weighted,support)
        y = acb_mat([[local[j,h]-sizes[j]*source[j,h] for h in (0,1)] for j in range(7)])
        gram = adj(local)*local
        delta = gram-adj(y)*y
        error = sum((2*norm(row)+d)*d for row,d in zip(y.tolist(),caps,strict=True))
        entries = [abs(z).lower() for z in delta.entries()]
        assert all(z.is_finite() for z in entries)
        value = (max(entries)-error).lower()
        assert value.is_finite() and value > arb(".14"), (support,value)
        mixed_gap = min(mixed_gap,value)
        # 正交投影不增大误差的 Frobenius 范数。投影后的 Gram 必须不大于局部 Gram。
        projected = y-projector*y
        error = (2*norm(projected.entries())+first_error)*first_error
        value = (-smaller_eigenvalue(gram-adj(projected)*projected)-error).lower()
        assert value.is_finite() and value > arb(".42"), (support,value)
        first_gap = min(first_gap,value)
        count += 1
    assert count == 729 and mixed_gap > arb(".14") and first_gap > arb(".42"), (count,mixed_gap,first_gap)
    print("每轨一站的支撑数、混合选线余量、全第一线余量",count,mixed_gap.str(12),first_gap.str(12),flush=True)


if __name__ == "__main__":
    root = Path(__file__).resolve().parent
    source = run_path(str(root/"ch39-fifth-full-positive-joint.py"))
    experiment(root,*source["verify"](root,return_operator="spectral"))
