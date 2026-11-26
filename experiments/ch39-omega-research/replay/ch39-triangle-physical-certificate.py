"""由三项原壳实读数的有限精度码，在整环面上认证五荷、两圈与曲率。"""
from itertools import product
from pathlib import Path
from runpy import run_path
import sys

from flint import acb, acb_mat, arb, ctx, fmpz_mat


root = Path(sys.argv[1]) if len(sys.argv) > 1 else Path(__file__).resolve().parent
physical = run_path(str(root / "ch39-first-three-mode-observables.py"))
joint = physical["joint"]
chiral = joint["chiral"]
geometry = chiral["source"]["geometry"]
faces = ((1, 6, 17), (1, 8, 9), (0, 7, 8))
bits = 26
codes = {8236: (65603085, 64036911, -66822133),
         8359: (67063892, 66230147, -65418989)}
delta_max = arb(23)/25*arb(2)**-20/14**2


def coefficients():
    edges = {(i, j): (distance, value) for i, j, distance, value in geometry["edges"]}
    result = []
    for (i, j, k), frequency in zip(faces, ((1, 0), (0, 1), (1, -1)), strict=True):
        a, b, c = (edges[e] for e in ((i, j), (j, k), (i, k)))
        observed = tuple(int(a[0] == d)+int(b[0] == d)-int(c[0] == d) for d in (724, 804))
        value = a[1]*b[1]*c[1].conjugate()
        if observed != frequency:
            assert tuple(-v for v in observed) == frequency
            value = value.conjugate()
        result.append(value)
    return result


def provenance(number, code, coeff):
    clocks = geometry["phases"](acb.zeta_zero(number).imag)
    phases = (clocks[724], clocks[804], clocks[724]/clocks[804])
    for k, z, p in zip(code, coeff, phases, strict=True):
        assert ((z*p).real*2**bits).floor().unique_fmpz() == k


def phase_box(code, coeff):
    # 整个量化格的两侧各放宽一个格宽；覆盖以实际值为中心的 ±2^-bits。
    values = [arb(k-1).union(arb(k+2))/2**bits for k in code]
    intervals = []
    for direction, target in zip(coeff[:2], values[:2], strict=True):
        angle = target.acos()
        offset = direction.arg()
        intervals.append((-offset-angle, -offset+angle))
    accepted, missed = [], []
    for x, y in product(*intervals):
        predicted = (coeff[2]*acb(0, x-y).exp()).real
        if predicted.overlaps(values[2]):
            accepted.append((x, y))
        else:
            missed.append(abs(predicted-values[2]).lower())
    assert len(accepted) == 1 and len(missed) == 3
    x, y = accepted[0]
    return x, y, min(missed)


def mass_bands(model, x, y):
    n = 19
    assert sum(1 for *_ij, d, _u in geometry["edges"] if d in (724, 804)) == 12
    incident = [sum(i == v or j == v for i, j, d, _u in geometry["edges"] if d in (724, 804))
                for v in range(n)]
    assert max(incident) == 4
    error = max(x.rad().upper(), y.rad().upper())
    phases = {724: acb(0, x.mid()).exp(), 804: acb(0, y.mid()).exp()}
    center = (15-delta_max/3).union(15)
    bands = {}
    for q in (0, *model["charges"]):
        matrix = acb_mat(n, n)
        for i, j, d, value in geometry["edges"]:
            matrix[i, i] += 1
            matrix[j, j] += 1
            matrix[i, j] = -(value*phases.get(d, 1))**q
            matrix[j, i] = matrix[i, j].conjugate()
        if q == 0:
            integer = fmpz_mat([[matrix[i, j].real.unique_fmpz() for j in range(n)] for i in range(n)])
            spectrum = [z for z, mult in integer.charpoly().complex_roots() for _ in range(mult)]
        else:
            spectrum = matrix.eig(multiple=True)
        assert all(z.imag.contains(0) for z in spectrum)
        levels = chiral["ordered"]([z.real for z in spectrum])
        if q == 0:
            assert levels[0].contains(0) and levels[1] > 1
            levels = levels[1:]
        # 实辐边的半正定减量至多 4δ，格内相位的厄米行和界为 4q·error。
        shift = arb(0, 4*q*error)-arb(0).union(4*delta_max)
        levels = [v+shift for v in levels]
        bands[q] = chiral["ordered"]([v*abs(v-center) for v in levels])
    return bands


def physical_windows(model, initial, bands):
    lo = max(bands[q][2].upper() for q in model["charges"])
    hi = min(bands[q][3].lower() for q in model["charges"])
    least_cut = 9+(1-delta_max)**2
    assert lo < least_cut and 10 < hi
    mu, _, alpha = joint["gauge"]["crossing_enclosure"](
        model, initial, bands, False, 724, "29.04", ".1184", recenter=True)
    assert arb(".117") < alpha < arb(".119")
    unit = 32*arb("2.435e18")/724**2
    even = [*bands[0], *(v for q in (2, 4, 6) for v in bands[q][3:])]
    decoupling = min((unit*h/mu).lower() for h in even)
    assert decoupling > 1
    vacuum = joint["vacuum"]
    squared = (arb("246.22")**2/(unit*bands[2][15]*arb("1.22089e19")))**2
    paired = vacuum["paired_spectrum"](vacuum["minimal"]["fields"](bands), min_gap="1e-18")
    rows = paired+vacuum["low_sector"](4, squared)
    left, right = (vacuum["newton"](rows, arb(x)) for x in (".03", ".033"))
    assert left < 1 < right
    coefficient = (9*arb(3).log()-12*arb(2).log())/2
    curve = squared**2*coefficient/(arb.pi()*arb("1e-122"))*(1-arb(0).union(squared/(arb(".03")*coefficient)))
    assert arb(".998") < curve < arb("1.002")
    return lo, hi, alpha, curve, decoupling, least_cut


def verify():
    model, initial = physical["model_inputs"]()
    ctx.prec = 192
    coeff = coefficients()
    assert acb.zeta_zero(1).imag > 14
    for number, code in codes.items():
        provenance(number, code, coeff)
        x, y, margin = phase_box(code, coeff)
        print("编号、角最大半宽与其余分支的第三读数间隔", number, x.rad(), y.rad(), margin)
        bands = mass_bands(model, x, y)
        lo, hi, alpha, curve, decoupling, least_cut = physical_windows(model, initial, bands)
        print("第三质量上界、第四质量下界、最轻偶重阈值除以交点", lo.str(12), hi.str(12), decoupling.str(10))
        print("强耦合、曲率上下端", *(w.str(12) for v in (alpha, curve) for w in (v.lower(), v.upper())))
        light_margin = min((least_cut-lo).lower(), hi-10)
        gauge_margin = min((alpha-arb(".117")).lower(), (arb(".119")-alpha).lower())
        curve_margin = min((curve-arb(".998")).lower(), (arb("1.002")-curve).lower())
        if number == 8236:
            assert bands[6][1] < 9 < bands[6][2]
            assert x.rad() < arb("1.066e-7") and y.rad() < arb("7.522e-8")
            assert margin > arb(".0503350")
            limits = (".58544", ".00068791", ".00106066", "1.4253")
        else:
            assert lo < 9 < hi
            assert x.rad() < arb("6.108e-7") and y.rad() < arb("1.389e-7")
            assert margin > arb(".0189259")
            limits = (".56754", ".00057819", ".00134511", "1.4294")
        assert all(read > arb(limit) for read, limit in
                   zip((light_margin, gauge_margin, curve_margin, decoupling), limits, strict=True))
        print("三窗离各边界的最小余量", *(v.str(11) for v in (light_margin, gauge_margin, curve_margin)))


if __name__ == "__main__":
    verify()
