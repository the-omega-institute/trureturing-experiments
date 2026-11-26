"""固定两个整数标尺，在联合零点检验球面单轨道变形引起的 Wilson 与谱矩变化。"""
from fractions import Fraction
from itertools import combinations, permutations
from pathlib import Path
from runpy import run_path
import sys

from flint import acb, acb_mat, arb, ctx


POLE = (17, 17, 17)
ORBIT = set(permutations((7, 17, 23)))
NUMBER = 2_306_888_100_085


def deform(points):
    ratio = Fraction(1001, 1000)
    h = (1-ratio**2)/(1+ratio**2)
    transverse = 2*ratio/(1+ratio**2)
    assert abs(h) < Fraction(1, 500) and transverse**2 == 1-h**2
    result = []
    for p in points:
        if p in ORBIT:
            z = Fraction(sum(p), 51)
            new = tuple((17*(z+h)+transverse*(coord-17*z))/(1+h*z) for coord in p)
        else:
            new = p
        assert sum(x*x for x in new) == 867
        result.append(new)
    assert sum(a != b for a, b in zip(result, points)) == 6
    return result


def enclose(v):
    return arb(v.numerator)/v.denominator


def transport(graph, x, y):
    dot, det = graph["dot"], graph["det"]
    re = 867*(867+dot(POLE, x)+dot(x, y)+dot(y, POLE))
    im = -det(POLE, x, y)
    z = acb(enclose(re), arb(867).sqrt()*enclose(im))
    return z/abs(z)


def verify(root):
    graph = run_path(str(root/"ch39-arc-relative-determinant.py"))
    shell = run_path(str(root/"ch39-three-ray-chiral-boundary.py"))
    read = run_path(str(root/"ch39-scale-free-geometric-pair.py"))["readout"]
    ctx.prec = 224
    points = graph["sites"]
    changed = deform(points)
    before = read(points)
    after = read(changed)
    assert before == after and before[0] == (32,724) and before[2] == (32,724,804)
    assert before[1] == {(i,j) for i,j,_,_ in graph["edges"]}
    height = acb.zeta_zero(NUMBER).imag
    phases = graph["phases"](height)
    pole, softened = shell["soft"]["pole"], shell["soft"]["star"]["orbit"]
    deficit = arb(23)/25*arb(2)**-20/(height**2+shell["soft"]["point"]["study"]["offset"])
    original = shell["laplacian"](1,height,phases=phases)
    matrix = acb_mat(19,19)
    edges = {}
    for i,j,d,native in graph["edges"]:
        stored = -1 if d == 200 else 1
        old = stored*transport(graph,points[i],points[j])
        new = stored*transport(graph,changed[i],changed[j])
        assert (old-native).contains(0)
        if d in (724,804):
            assert changed[i] == points[i] and changed[j] == points[j]
        clock = phases.get(d,acb(1))
        spoke = pole in (i,j)
        if spoke:
            assert (new-1).contains(0) and (clock-1).contains(0)
        w = 1-deficit if spoke and (i in softened or j in softened) else arb(1)
        old_phase, new_phase = old*clock, new*clock
        edges[i,j] = (old_phase, new_phase, w)
        matrix[i,i] += w
        matrix[j,j] += w
        matrix[i,j] = -w*new_phase
        matrix[j,i] = -w*new_phase.conjugate()
    def eigenvalues(m):
        values = m.eig(multiple=True)
        assert len(values) == 19 and all(z.imag.contains(0) for z in values)
        levels = sorted((z.real for z in values), key=float)
        assert all(a < b for a, b in zip(levels, levels[1:]))
        return levels

    before_levels, after_levels = eigenvalues(original), eigenvalues(matrix)
    assert arb("1.3696573220e-5") < after_levels[0]-before_levels[0] < arb("1.3696573221e-5")
    print("最低荷一特征值与位移", before_levels[0].str(17), after_levels[0].str(17), (after_levels[0]-before_levels[0]).str(17))
    shifts = acb(0)
    witness = tuple(points.index(p) for p in ((1,29,5), (23,13,13), (7,17,23)))
    assert witness == (1,8,9)
    count = 0
    for i,j,k in combinations(range(19),3):
        if not all(e in edges for e in ((i,j),(j,k),(i,k))):
            continue
        a,b,c = (edges[e] for e in ((i,j),(j,k),(i,k)))
        values = [a[n]*b[n]*c[n].conjugate() for n in (0,1)]
        shifts += a[2]*b[2]*c[2]*(values[1]-values[0]).real
        if (i,j,k) == witness:
            shift = (values[1]-values[0]).real
            assert arb("-2.5229723732e-5") < shift < arb("-2.5229723730e-5")
            print("指定804钟三角面原、新实部及位移",values[0].real.str(18),values[1].real.str(18),shift.str(18))
        count += 1
    assert count == 529 and all(e in edges for e in combinations(witness,2))
    change = matrix**3-original**3
    trace = sum(change[i,i] for i in range(19))
    assert trace.imag.contains(0) and (trace+6*shifts).contains(0)
    assert arb(".00561804236851") < trace.real < arb(".00561804236852")
    print("矩阵三阶迹变化、三角和重算", trace.real.str(18),(-6*shifts.real).str(18))
    print("轨道代表变形后",changed[points.index((7,17,23))])


if __name__ == "__main__":
    verify(Path(sys.argv[1]) if len(sys.argv)>1 else Path(__file__).resolve().parent)
