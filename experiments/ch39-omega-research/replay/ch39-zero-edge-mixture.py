"""时钟 1 的带荷边读数与全矩阵逐项误差下的最优全谱混合。"""
from itertools import combinations, permutations
from pathlib import Path
from runpy import run_path

from flint import acb, acb_mat, arb, arb_poly, ctx

memory = run_path(str(Path(__file__).with_name("ch39-three-clock-minimal-memory.py")))


def phase_tolerance(left, target, right):
    chord = left-right
    normal = acb(0, 1)*chord/abs(chord)
    if (normal.conjugate()*(target-right)).real < 0:
        normal = -normal
    rotated = normal.conjugate()*target
    m, k = rotated.real, rotated.imag
    boundary = (normal.conjugate()*right).real
    distance = m-boundary
    assert m > 0 and boundary > 0 and distance > 0 and k < 0
    # u=tan(|角误差|/2)。有理化的正根避免两近数相减。
    u = distance/(abs(k)+(k*k+(m+boundary)*distance).sqrt())
    assert 0 < u < 1
    box = arb(0).union(-u)
    def factor(t):
        return acb(1, t)/acb(1, -t)
    projection = ((factor(box)*target-right)/chord).real
    assert 0 < projection < 1
    endpoint = (factor(-u)*target-right)/chord
    assert endpoint.imag.contains(0) and 0 < endpoint.real < 1
    angle = 2*u.atan()
    assert arb("3.09312498382782e-12") < angle < arb("3.09312498382783e-12")
    print("First relative phase ambiguity (rad):", angle.str(24),
          "left-neighbor weight:", endpoint.real.str(24))


def verify():
    ctx.prec = 384
    _, _, _, g0, g1, g2 = memory["coefficients"]()
    model = memory["REFINE"]
    nodes, edges = model["geometry"]()
    base = acb_mat(memory["identity"](19)-g0-g1-g2)
    charged = []
    for i, j, weight, exponent in edges:
        if j >= 19:
            continue
        weight = model["ball"](weight)
        if exponent:
            charged.append((i, j, weight, exponent))
        else:
            base[i, j] -= weight
            base[j, i] -= weight
    assert len(charged) == 6

    def pencil(s):
        unit = acb(1, s)
        w = unit/unit.conjugate()
        result = acb_mat(base)
        for i, j, weight, charge in charged:
            result[i, j] -= weight*w**charge
            result[j, i] -= weight*w**(-charge)
        return result

    edge = (4, 6)
    ts = [arb(i)/8 for i in range(-6, 7)]
    polys = [arb_poly() for _ in range(3)]
    for i, t in enumerate(ts):
        matrix = pencil(t)
        den = matrix.det().real*(1+t*t)**6
        num = matrix.inv()[edge]*den
        basis = arb_poly([1])
        for j, other in enumerate(ts):
            if i != j:
                basis *= arb_poly([-other, 1])*(1/(ts[i]-other))
        for k, value in enumerate((den, num.real, num.imag)):
            polys[k] += basis*value
    q, a, b = polys
    assert acb.zeta_zero(1).imag > 10
    domain = arb(0).union(arb(1)/20)
    sa, sb = [p.derivative()*q-p*q.derivative() for p in (a, b)]
    assert q(domain) > 0 and 0 not in sb(domain)
    curvature = sa*sb.derivative()-sb*sa.derivative()
    assert 0 not in curvature(domain)
    # 反向运行与正向运行共用端口参考。所有混合条目均落在此凸锥内。
    assert a(domain) > 2*abs(b(domain))
    assert 2*(arb(1)/2).atan()+arb("4e-12") < arb.pi()/2

    left, target, right = [pencil(1/(2*acb.zeta_zero(n).imag)).inv()
                           for n in (11430, 11431, 11432)]
    chord, offset = left[edge]-right[edge], right[edge]-target[edge]
    weight = -(chord.conjugate()*offset).real/abs(chord)**2
    gap = abs((chord.conjugate()*offset).imag)/abs(chord)
    assert 0 < weight < 1 and gap > 0

    directed = {(a, b): (w, q) for i, j, w, e in edges
                for a, b, q in ((i, j, e), (j, i, -e))}
    orbit = set()
    for perm in permutations(range(3)):
        act = [nodes.index(tuple(p[k] for k in perm)) for p in nodes]
        sign = (-1)**sum(perm[i] > perm[j] for i, j in combinations(range(3), 2))
        assert all(directed[act[i], act[j]] == (w, sign*e) for (i, j), (w, e) in directed.items())
        orbit.add(tuple(sorted((act[edge[0]], act[edge[1]]))))
    assert len(orbit) == 6
    difference = right-target+weight*(left-right)
    for i in range(19):
        for j in range(i, 19):
            error = abs(difference[i, j])
            if (i, j) in orbit:
                assert error.overlaps(gap)
            else:
                assert error < gap
    phase_tolerance(left[edge], target[edge], right[edge])
    print("Optimal maximum entry gap at clock 1:", gap.str(25),
          "weight on 11430:", weight.str(25))


if __name__ == "__main__":
    verify()
