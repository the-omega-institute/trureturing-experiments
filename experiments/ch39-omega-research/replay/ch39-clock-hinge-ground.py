"""实际零点弦钟给出的两类面积、四／八次计数的占据能隙与阈值。"""
from math import isqrt
from pathlib import Path
from runpy import run_path

from flint import acb, arb, ctx


def matched_geometry(n, k):
    """两个计数分别是 k 与 n-k。"""
    other = n - k
    ratio = arb(k) / other
    return geometry_at_ratio(n, other, ratio)


def geometry_at_ratio(n, k6, ratio):
    t = 2*(1-ratio**2) + (4*(1-ratio**2)**2 + ratio**2).sqrt()
    b = 4*k6 / (n*(4*t-1).sqrt())
    four = 2*arb.pi() - 4*((3*t-2)/(2*(1+t))).acos()
    six = 2*arb.pi() - 6*((3-2*t)/(2*(1+t))).acos()
    volume = arb(5).sqrt()*b**2*(3*t-t**2-1)/4
    curvature = 2*arb(k6)/n*(30*ratio*four + 20*six)/volume
    assert curvature < 0
    return t, four, six, curvature


def critical_strength(n, k, deficit, m):
    assert 3 <= k < n and (deficit > 0 or deficit < 0)
    jump = 2*k - 1 if deficit > 0 else 2*k + 1
    return abs(deficit)*n**(4*m) / (4*arb.pi()*jump**(2*m))


def largest_unit_strength_size(probabilities, target):
    n0 = 5432
    epsilon = arb("0.0000001")
    # 单位强度的必要角亏界随 n 递减。反解前用严格单调性夹住 p。
    limit = 4*arb.pi()*(2*n0+1)**2/n0**4
    for p in (target-epsilon, target+epsilon):
        assert abs(geometry_at_ratio(1, 1-p, p/(1-p))[1]) > limit
    radius = epsilon + arb(1)/(2*n0)
    candidates = []
    for index, p in enumerate(probabilities, 10855):
        if abs(p-target) < radius:
            candidates.append(index)
        else:
            assert abs(p-target) > radius
    assert candidates == [11366]
    p = probabilities[11366 - 10855]
    assert abs(p-target) > epsilon + arb(1)/(2*5800)
    hits = []
    for n in range(n0, 5800):
        k = int((n*p+arb(1)/2).floor().unique_fmpz())
        _, four, six, curvature = matched_geometry(n, k)
        strength = max(critical_strength(n, k, four, 1), critical_strength(n, n-k, six, 1))
        if strength < 1:
            hits.append((n, k))
            assert -arb("5.811e-14") < curvature < -arb("5.810e-14")
            print("单位四读强度可达的最大槽位、零点、负载、临界强度", n, 11366, k,
                  strength.str(12))
        else:
            assert strength > 1
    assert hits == [(5432, 2911)]


def clock_window():
    base = run_path(str(Path(__file__).with_name("ch39-arc-relative-determinant.py")))
    ctx.prec = 192
    zeros = acb.zeta_zeros(10855, 1238)
    assert len(zeros) == 1238 and all(a.imag < b.imag for a, b in zip(zeros, zeros[1:]))
    probabilities = [(1 + base["phases"](z.imag)[724].real)/2 for z in zeros]
    lower = (3-arb(5).sqrt())/2
    upper = 1-lower
    assert arb(".16936") < probabilities[11431-10855] < arb(".16937")
    target = 4 - 2*arb(3).sqrt()
    assert all(abs(p-target) > arb(".00008655") for p in probabilities)
    # N≥10000 时舍入误差最多 0.00005，两个端点给四价角亏的保守下界。
    for sign in (-1, 1):
        p = target + sign*arb(".00003655")
        _, four, _, _ = geometry_at_ratio(1, 1-p, p/(1-p))
        assert abs(four) > arb(".00087")
    largest_unit_strength_size(probabilities, target)
    for n, count, expected_load in [(1000, 184, 536), (2**20, 183, 562021)]:
        rows = []
        for index, probability in enumerate(probabilities, 10855):
            rounded = (n*probability + arb(1)/2).floor().unique_fmpz()
            assert rounded is not None
            k = int(rounded)
            p = arb(k)/n
            if p < lower or p > upper:
                continue
            assert lower < p < upper
            t, four, six, curvature = matched_geometry(n, k)
            strength = [max(critical_strength(n, k, four, m),
                            critical_strength(n, n-k, six, m)) for m in (1, 2)]
            rows.append((index, k, t, abs(curvature), *strength))
        best = min(rows, key=lambda row: float(row[3]))
        ties = [row[0] for row in rows if row[1] == best[1]]
        assert all(row[1] == best[1] or row[3] > best[3] for row in rows)
        assert (len(rows), ties, best[1]) == (count, [11366], expected_load)
        assert all(row[1] == best[1] or (row[4] > best[4] and row[5] > best[5])
                   for row in rows)
        print("N、可实现的编号数、最小曲率绝对值的编号、四价负载", n, len(rows), ties, best[1])
        print("形状、曲率绝对值、四读和八读临界强度", *[v.str(12) for v in best[2:]])


def pell_limit():
    ctx.prec = 512
    w, k4 = 2, 1
    for _ in range(12):
        w, k4 = 7*w + 12*k4, 4*w + 7*k4
    k6 = w//2
    assert 3*k4**2 - 4*k6**2 == -1
    n = 1 + isqrt(20*k6**2//3)
    _, four, six, _ = geometry_at_ratio(n, k6, arb(k4)/k6)
    assert four > 0 and six < 0
    for m in (1, 2):
        threshold = critical_strength(n, k6, six, m)
        assert threshold > critical_strength(n, k4, four, m)
        limit = 5*arb(3).sqrt() / (4*arb.pi()*(arb(3)/5)**m)
        assert abs(threshold/n**(2*m-2)-limit) < arb("1e-11")
        print("读数次数与 Pell 支归一临界强度极限", 4*m, limit.str(13))
    n = k4 + k6
    q = 2*arb(3).sqrt() - 3
    _, four, six, curvature = matched_geometry(n, k4)
    curvature_limit = 27*(arb(5)/3).sqrt()/(32*q**5)
    assert abs(abs(curvature)*n**4-curvature_limit) < arb("1e-10")
    for m in (1, 2):
        limit = 3*arb(3).sqrt()/(16*arb.pi()*q**2*(2*q)**(2*m))
        assert abs(critical_strength(n, k6, six, m)/n**(2*m-2)-limit) < arb("1e-11")
        print("互补计数下全谱可实现的临界强度下极限", 4*m, limit.str(13))


if __name__ == "__main__":
    clock_window()
    pell_limit()
