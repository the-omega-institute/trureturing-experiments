"""十九站原边表上，核验共同钟移的两个极小极大角误差及其达到构造。"""
from itertools import permutations
from pathlib import Path
from runpy import run_path

from flint import acb, arb


def angle(t):
    return acb(0, t).exp().arg()


def verify():
    shell = run_path(str(Path(__file__).with_name("ch39-arc-relative-determinant.py")))
    edges = [(i,j,d) for i,j,d,_ in shell["edges"]]
    outer = set(range(6))
    inner = {6, 7, 8}
    endpoints = {i: {d: j for k,j,d in edges if k == i and d in (724,804)}
                 for i in outer}
    assert {(v[724],v[804]) for v in endpoints.values()} == set(permutations(inner,2))
    assert {(i,j) for i,j,d in edges if d==544} == {(i,14) for i in outer}
    for i,j,d in edges:
        if d>372:
            assert int(j in inner)-int(i in inner) == int(d in (724,804))
    tails = [(i,j) for i,j,d in edges if 544<d<724]
    assert len(tails) == len({j for _,j in tails}) == 6
    assert all(i in outer and j not in outer|inner for i,j in tails)

    heights = [acb.zeta_zero(n).imag for n in (2345666,964890)]
    delta = abs(angle((shell["theta"](724)-shell["theta"](804))*(heights[0]-heights[1])))
    pi = arb.pi()
    assert 2*pi/3 < delta < pi
    r0, r1 = delta/2, 2*pi/3-delta/2
    assert arb("1.4660591778386") < r0 < arb("1.4660591778387")
    assert arb(".62833592455456") < r1 < arb(".62833592455459")
    constant = [r0 if j in inner else arb(0) for j in range(19)]
    triangular = [arb(0) for _ in range(19)]
    for k,j in enumerate(sorted(inner)):
        triangular[j] = 2*pi*k/3
    for i, ends in endpoints.items():
        left, right = triangular[ends[724]]-delta, triangular[ends[804]]
        triangular[i] = left + angle(right-left)/2
    for i,j in tails:
        triangular[j] = triangular[i]

    for cutoff, phases, radius in [(372, constant, r0), (544, triangular, r1)]:
        errors = [abs(angle((delta if d==724 else 0)+phases[i]-phases[j]))
                  for i,j,d in edges if d>cutoff]
        assert all(e < radius + arb("1e-40") for e in errors)
        assert any(e > radius - arb("1e-40") for e in errors)
        print("平方弦长上限、达到的角误差",cutoff,radius.str(16))




def verify_department_charges():
    shell = run_path(str(Path(__file__).with_name("ch39-arc-relative-determinant.py")))
    up, down = [acb.zeta_zero(n).imag for n in (2345666,964890)]
    inner = (6,7,8)
    pi = arb.pi()
    targets = {}
    clock_sites = {i:{} for i in range(6)}
    tails = []
    for i,j,d,z in shell["edges"]:
        if d<=544:
            continue
        value = z**-6
        if d in (724,804):
            value *= acb(0,(shell["theta"](d)-shell["theta"](244))*(-4*up-2*down)).exp()
            clock_sites[i][d] = j
        else:
            assert d==632
            tails.append((i,j))
        targets[i,j] = value.arg()
    assert len(tails) == len({j for _,j in tails}) == 6
    # 三个内点按 6→7→8→6 定向，三个角差都可取 2π/3。
    phases = [arb(0) for _ in range(19)]
    for k,i in enumerate(inner):
        phases[i] = -2*pi*k/3
    path_angles = {}
    for i,ends in clock_sites.items():
        j, k = ends[724],ends[804]
        target = angle(targets[i,j]-targets[i,k])
        path_angles[j,k] = target
        left, right = phases[j]-targets[i,j], phases[k]-targets[i,k]
        phases[i] = left+angle(right-left)/2
    for i,j in tails:
        phases[j] = phases[i]+targets[i,j]
    cycle = [(6,7),(7,8),(8,6)]
    # 三次路径绕过六条边，任意规范下最大边误差至少是环量角距的六分之一。
    lower = abs(angle(sum(path_angles[pair] for pair in cycle)))/6
    upper = arb(0)
    for (i,j),target in targets.items():
        upper = upper.max(abs(angle(target+phases[i]-phases[j])))
    assert arb(".39645665403697162") < lower
    assert upper < arb(".39645665403697164")
    print("两部门荷数和背景同时保留时的环量下界、构造上界",lower.str(18),upper.str(18))


if __name__ == "__main__":
    verify()
    verify_department_charges()
