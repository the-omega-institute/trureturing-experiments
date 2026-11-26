"""非实时钟的单个复对角读数区分 11431 与全正零点谱上的非相干混合。"""
from itertools import combinations, permutations
from pathlib import Path
from runpy import run_path

from flint import acb, acb_mat, acb_poly, arb, arb_poly, ctx

single = run_path(str(Path(__file__).with_name("ch39-zero-single-port.py")))
memory = single["memory"]


def interpolate():
    _, _, _, g0, g1, g2 = memory["coefficients"]()
    nodes, edges = memory["REFINE"]["geometry"]()
    directed = {(a, b): (w, q) for i, j, w, exponent in edges
                for a, b, q in ((i, j, exponent), (j, i, -exponent))}
    for perm in permutations(range(3)):
        action = [nodes.index(tuple(point[h] for h in perm)) for point in nodes]
        sign = (-1)**sum(perm[i] > perm[j] for i, j in combinations(range(3), 2))
        assert all(directed[action[i], action[j]] == (w, sign*q)
                   for (i, j), (w, q) in directed.items())
    shell = memory["REFINE"]
    old = [edge for edge in edges if edge[1] < 19]
    z = acb(-1, arb(3).sqrt())/2
    base = z*acb_mat(memory["identity"](19))-z**2*acb_mat(g0)-z*acb_mat(g1)-acb_mat(g2)
    xs, rows = [], []
    for n in range(1, 8):
        w = arb(n)
        xs.append((w+1/w)/2-1)
        pencil = acb_mat(base)
        # 插值点取正实 w，反向带荷因子为 w 的倒数。
        for i, j, weight, exponent in old:
            value = shell["ball"](weight)
            pencil[i, j] -= value*w**exponent
            pencil[j, i] -= value*w**(-exponent)
        det = pencil.det()
        inv = pencil.inv()
        rows.append([det, *(det*inv[i, i] for i in range(19))])
    result = [acb_poly() for _ in rows[0]]
    for i, row in enumerate(rows):
        basis = arb_poly([1])
        for j, other in enumerate(xs):
            if i != j:
                basis *= arb_poly([-other, 1])*(1/(xs[i]-other))
        for k, value in enumerate(row):
            result[k] += acb_poly(basis)*value
    return result


def parts(poly):
    return (arb_poly([c.real for c in poly.coeffs()]),
            arb_poly([c.imag for c in poly.coeffs()]))


def verify():
    ctx.prec = 440
    den, *nums = interpolate()
    qr, qi = parts(den)
    q = qr*qr + qi*qi
    domain = arb(-1).union(arb(0))/200
    assert acb.zeta_zero(1).imag > 10 and q(domain) > 0
    xs = [-1/(2*(acb.zeta_zero(n).imag**2+arb(1)/4)) for n in (11430,11431,11432)]
    gaps, weights = [], []
    for num in nums:
        pr, pi = parts(num)
        a, b = pr*qr+pi*qi, pi*qr-pr*qi
        sa, sb = [p.derivative()*q-p*q.derivative() for p in (a, b)]
        curvature = sa*sb.derivative()-sb*sa.derivative()
        assert 0 not in sa(domain) and 0 not in curvature(domain)
        left, target, right = [num(x)/den(x) for x in xs]
        chord, offset = left-right, right-target
        weight = -(chord.conjugate()*offset).real/abs(chord)**2
        assert 0 < weight < 1
        gap = abs((chord.conjugate()*offset).imag)/abs(chord)
        assert gap > 0 and gap.overlaps(abs(offset+weight*chord))
        gaps.append(gap)
        weights.append(weight)
    best = max(range(19), key=lambda i: float(gaps[i]))
    tied = [i for i in range(19) if gaps[i].overlaps(gaps[best])]
    points = [(13,13,23), (13,23,13), (23,13,13)]
    sites = sorted(memory["REFINE"]["LOCAL"]["shell_grid"](1) | {(17,17,17)})
    assert tied == [sites.index(point) for point in points]
    assert all(gaps[i] < gaps[best] for i in range(19) if i not in tied)
    print("Largest single-port separation, sites", tied,
          "gap", gaps[best].str(22), "left weight", weights[best].str(22))


if __name__ == "__main__":
    verify()
