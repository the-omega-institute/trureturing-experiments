"""用十五个固定参考方向包围整个双钟环面的原壳轻线，核验粗网格的投影链界。"""
from pathlib import Path
from runpy import run_path
import sys

from flint import acb, acb_mat, arb, ctx


def adjoint(matrix):
    return matrix.conjugate().transpose()


def rotation_bases(shell):
    sites = shell["geo"]["sites"]
    rotation = [sites.index((y, z, x)) for x, y, z in sites]
    assert rotation[shell["soft"]["pole"]] == shell["soft"]["pole"]
    tuned = set(shell["soft"]["star"]["orbit"])
    assert len(tuned) == 3 and {rotation[i] for i in tuned} == tuned
    omega = acb(-1, arb(3).sqrt())/2
    for k in range(3):
        vectors = []
        for i in range(19):
            orbit = (i, rotation[i], rotation[rotation[i]])
            if min(orbit) != i:
                continue
            v = [acb(0)]*19
            if rotation[i] == i:
                if k:
                    continue
                v[i] = 1
            else:
                for j, m in enumerate(orbit):
                    v[m] = omega**(-k*j)/arb(3).sqrt()
            vectors.append(v)
        yield acb_mat(vectors).transpose()


def clock_matrices(shell, q):
    base = acb_mat(19,19)
    links = [acb_mat(19,19), acb_mat(19,19)]
    ends = [[], []]
    pole = shell["soft"]["pole"]
    sites = shell["geo"]["sites"]
    rotation = [sites.index((y,z,x)) for x,y,z in sites]
    for i,j,d,z in shell["geo"]["edges"]:
        base[i,i] += 1
        base[j,j] += 1
        entry = -acb(1) if pole in (i,j) else -z**q
        if d in (724,804):
            which = int(d == 804)
            assert rotation[i] < rotation[j]
            links[which][i,j] = entry
            ends[which].append((i,j))
        else:
            base[i,j] = entry
            base[j,i] = entry.conjugate()
    for edges in ends:
        sources, targets = zip(*edges)
        assert len(set(sources)) == len(edges) == 6
        assert len(set(targets)) == 3 and all(targets.count(t)==2 for t in targets)
    assert set(i for edges in ends for i,j in edges).isdisjoint(j for edges in ends for i,j in edges)
    return base, links


def first_line(matrix):
    eig, vectors = matrix.eig(right=True, algorithm="rump")
    order = sorted(range(len(eig)), key=lambda j:float(eig[j].real))
    lo = eig[order[0]].real
    assert all(z.imag.contains(0) for z in eig) and all(lo < eig[j].real for j in order[1:])
    gap = min((eig[j].real-lo).lower() for j in order[1:])
    v = [vectors[i,order[0]] for i in range(matrix.nrows())]
    norm = sum(abs(z)**2 for z in v).sqrt()
    return [z/norm for z in v], gap


def verify(root):
    shell = run_path(str(root/"ch39-three-ray-chiral-boundary.py"))
    ctx.prec = 144
    bases = list(rotation_bases(shell))
    grid = 32
    roots = [acb(0, 2*arb.pi()*i/grid).exp() for i in range(grid)]
    # 格心至全格的算子位移，另计三条辐边的独立上界。
    change = 4*arb(2).sqrt()*(arb.pi()/(2*grid)).sin() + 6*arb("1e-8")
    worst = arb(0)
    for q in (1,-4,2,-3,6):
        base, links = clock_matrices(shell, q)
        angles = []
        for Q in bases:
            d = adjoint(Q)*base*Q
            anchor, _ = first_line(d)
            first, second = [adjoint(Q)*m*Q for m in links]
            a_terms = [first*z + adjoint(first*z) for z in roots]
            b_terms = [second*z + adjoint(second*z) for z in roots]
            upper = arb(0)
            for part in a_terms:
                center = d+part
                for other in b_terms:
                    v, gap = first_line(center+other)
                    overlap = sum(z.conjugate()*w for z,w in zip(anchor,v))
                    sq = abs(overlap)**2
                    assert sq > 0 and gap > 1+2*change
                    midpoint_sine = (1-sq.lower()).max(arb(0)).sqrt()
                    limit = (midpoint_sine + change/(gap-change)).upper()
                    upper = max(upper, limit)
            assert upper < arb(1)/3
            angles.append(upper)
        worst = max(worst, *angles)
        print("荷、三个特征标全域的主角正弦上界",q,*[v.str(11) for v in angles],flush=True)
    bound = arb(1)/3
    probability = (1-2*bound**2)**2
    phase = (bound**2/(1-bound**2)).asin()
    assert 25*phase < arb.pi() and (probability-arb(49)/81).contains(0)
    width = arb(1)/16
    size = 2*arb(2).sqrt()*(width.exp()-1)
    resolvent = 1/(arb(1)/2-size)
    projection = resolvent/2
    gradient = resolvent**2*arb(2).sqrt()*width.exp()/2
    curvature = 2*projection*gradient**2
    assert resolvent < 4 and curvature < 800
    theta = shell["geo"]["theta"]
    ratio = (theta(804)-theta(244))/(theta(724)-theta(244))
    mismatch = abs(arb(941)/829-ratio)
    periodic = 1600/((arb(829+941)/16).exp()-1)
    assert mismatch < arb("5.52e-7") and periodic < arb("1.45e-45")
    print("有理斜率 941/829 的钟速比误差、有效钟角坐标下的整周期平均曲率界",
          mismatch.str(12), periodic.str(12))
    print("双复角带宽、Berry 曲率在复带的统一绝对值上界",width,curvature.str(12))
    print("全域最大证书上界、任意两背景的投影成功率下界、正相位规范下每边主辐角绝对值上界",
          worst.str(12), probability, phase.str(12))


if __name__ == "__main__":
    root = Path(sys.argv[1]).resolve() if len(sys.argv)>1 else Path(__file__).resolve().parent
    del sys.argv[1:]
    verify(root)
