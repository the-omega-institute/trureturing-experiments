"""Certify the winner fan of the parity-weighted vacuum family of the 38-mode shell.

For every zero n in the window 10815..12092 the 38-mode particle-hole shell gives
the many-body vacuum energy E(n), the lowest quasiparticle energy g(n) and the
Pfaffian sector s(n) in {+1,-1}.  The family studied is

    F(n; mu_plus, mu_minus) = E(n) + mu_{s(n)} g(n),        (mu_plus, mu_minus) in R^2.

(0,0) is the vacuum energy, the diagonal mu -> -infinity the parity gap, (0,1) and
(1,0) the two fixed-parity ground energies.  The script certifies the exact
lower envelopes, the set of attainable winners, the monotone no-go for 10875,
11737 and 11488, and the minimal unification cost of the two native criteria.
"""
import math
from pathlib import Path
from runpy import run_path

import numpy as np
from flint import acb, acb_mat, arb, ctx

root = Path(__file__).parent
model = run_path(str(root / "ch39-majorana-vacuum-energy.py"))
parity = run_path(str(root / "ch39-fixed-parity-equilibrium.py"))
shell = model["shell"]
unit = acb_mat([[int(i == j) for j in range(38)] for i in range(38)])
FIRST, COUNT = 10815, 1278


def frobenius(m):
    return sum(abs(z).upper()**2 for z in m.entries()).sqrt().upper()


def centre(m):
    return acb_mat([[m[i, j].mid() for j in range(m.ncols())] for i in range(m.nrows())])


def spectrum(h):
    """Certified E (half the sum of the 19 lowest eigenvalues) and g (20th eigenvalue)."""
    _, approx = np.linalg.eigh(np.array(h.tolist(), dtype=complex))
    v = acb_mat(approx.tolist())
    # Proposal only: orthonormalise, then one first-order Jacobi correction in
    # ball precision.  The certificate below does not depend on how v was found.
    v = centre(v * (3 * unit - v.conjugate().transpose() * v) / 2)
    a = v.conjugate().transpose() * h * v
    d = [a[i, i].real.mid() for i in range(38)]
    x = acb_mat([[0 if i == j else (a[i, j] / (d[j] - d[i])).mid() for j in range(38)]
                 for i in range(38)])
    v = centre(v * (unit + x))
    for _ in range(4):
        v = v * (3 * unit - v.conjugate().transpose() * v) / 2
        adj = v.conjugate().transpose()
        defect = frobenius(adj * v - unit)
        if defect < arb("1e-24"):
            break
        v = centre(v)
    assert defect < arb("1e-24")
    hv = h * v
    ray = adj * hv
    diag = [ray[i, i].real.mid() for i in range(38)]
    residual = acb_mat([[hv[i, j] - v[i, j] * diag[j] for j in range(38)] for i in range(38)])
    # Same perturbation bound as the vacuum certificate: ||H - Q D Q^*|| <= error,
    # Q the unitary polar factor of v.
    error = 2 * frobenius(residual) + defect * (
        frobenius(h) + 3 * max(abs(z).upper() for z in diag))
    assert error < arb("1e-12")
    d = sorted(diag)
    # Weyl: the k-th eigenvalue lies within error of d[k].  Particle-hole symmetry
    # pairs the spectrum; the 20th eigenvalue is certified positive.
    for k in range(19):
        assert abs(d[k] + d[37 - k]) < 2 * error
    g = d[19] + arb(0, error)
    assert g > 0
    energy = (sum(d[:19]) + arb(0, 19 * error)) / 2
    return energy, g


def collect():
    ctx.prec = 128
    data = {}
    for n, z in enumerate(acb.zeta_zeros(FIRST, COUNT), FIRST):
        parent = shell["parent_data"](z.imag)
        energy, g = spectrum(model["paired_hamiltonian"](*parent))
        data[n] = (energy, g, parity["vacuum_sector"](parent))
    return data


def orient(a, b, c):
    return (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0])


def envelope(points):
    """Certified classification of argmin_n y_n + c x_n over all real c.

    Returns the vertex chain v_0..v_k (x increasing) and the breakpoints
    c_0 > ... > c_{k-1}; v_0 wins on (c_0, inf), v_i on (c_i, c_{i-1}),
    v_k on (-inf, c_{k-1}), each uniquely; at c_i exactly v_i and v_{i+1} tie.
    """
    mids = sorted(points, key=lambda n: (float(points[n][0].mid()), float(points[n][1].mid())))
    chain = []
    for n in mids:
        while len(chain) >= 2 and float(orient(points[chain[-2]], points[chain[-1]], points[n]).mid()) <= 0:
            chain.pop()
        chain.append(n)
    first, last = chain[0], chain[-1]
    assert all(points[m][0] > points[first][0] for m in points if m != first)
    assert all(points[m][0] < points[last][0] for m in points if m != last)
    for a, b in zip(chain, chain[1:]):
        assert points[b][0] > points[a][0]
        # Every other point lies strictly above the supporting line of each edge;
        # with the next vertex this also certifies c_{i-1} > c_i.
        for m in points:
            if m not in (a, b):
                assert orient(points[a], points[b], points[m]) > 0, (a, b, m)
    breaks = [-(points[b][1] - points[a][1]) / (points[b][0] - points[a][0])
              for a, b in zip(chain, chain[1:])]
    return chain, breaks


def value(points, n, c):
    return points[n][1] + c * points[n][0]


def interval_bounds(values):
    lo = np.array([math.nextafter(float(v.lower()), -math.inf) for v in values])
    hi = np.array([math.nextafter(float(v.upper()), math.inf) for v in values])
    assert np.all(lo <= hi)
    return lo, hi


def verify():
    data = collect()
    ns = sorted(data)
    energy = {n: data[n][0] for n in ns}
    gap = {n: data[n][1] for n in ns}
    sector = {n: data[n][2] for n in ns}
    plus = [n for n in ns if sector[n] == 1]
    minus = [n for n in ns if sector[n] == -1]
    assert (len(plus), len(minus)) == (611, 667)
    # Consistency with the published vacuum and gap certificates.
    for n in (11618, 11431):
        old = model["vacuum"](model["hamiltonian"](acb.zeta_zero(n).imag))
        assert old.overlaps(energy[n]) and energy[n].rad() < old.rad()
    assert arb("0.0236718") < gap[11431] < arb("0.0236719")
    assert sector[11177] == 1 and sector[11618] == -1

    # --- one-parameter envelopes -------------------------------------------------
    diagonal = {n: (gap[n], energy[n]) for n in ns}
    chain, breaks = envelope(diagonal)
    assert chain == [11681, 10928, 12003, 10931, 11601, 11874, 11202, 11474, 11075,
                     11746, 11346, 11618, 11619, 11891, 11491, 11092, 11363, 11431]
    expected = ["494.0415", "0.94514", "0.93997", "0.90784", "0.77395", "0.71762", "0.60284",
                "0.45398", "0.44907", "0.29482", "0.10538", "-0.06733", "-0.11466",
                "-0.40180", "-1.03760", "-2.61783", "-7.42210"]
    for b, e in zip(breaks, expected, strict=True):
        assert abs(b - arb(e)) < arb("0.0006"), (b, e)
    # Named diagonal points: mu=0 vacuum energy, mu=1 opposite-parity level, mu=-1.
    named = {0: 11618, 1: 10928, -1: 11491}
    for mu, winner in named.items():
        assert all(value(diagonal, m, mu) > value(diagonal, winner, mu) for m in ns if m != winner)
    print("diagonal (parity-blind) winners, mu from +inf to -inf:", chain)
    print("  breakpoints", [b.str(6) for b in breaks])

    classes = {}
    for label, members in (("+", plus), ("-", minus)):
        points = {n: diagonal[n] for n in members}
        classes[label] = (points, *envelope(points))
    assert classes["+"][1] == [11681, 10928, 11321, 11993, 11050, 11448, 11177,
                               11576, 11304, 11977, 11033, 11431]
    assert classes["-"][1] == [11531, 11532, 12001, 11731, 12002, 12003, 10931, 11601,
                               11874, 11202, 11474, 11075, 11746, 11346, 11618, 11619,
                               11891, 11491, 11092, 11363, 10965]
    for label in "+-":
        print("sector", label, "winners:", classes[label][1])
        print("  breakpoints", [b.str(6) for b in classes[label][2]])

    # Every sector vertex is the unique global winner at an explicit (mu+, mu-).
    big = arb(10)**5
    attainable = []
    for label, other in (("+", "-"), ("-", "+")):
        points, vertices, cuts = classes[label]
        for i, v in enumerate(vertices):
            hi = cuts[i - 1].mid() if i > 0 else cuts[0].mid() * 2 + 1
            lo = cuts[i].mid() if i < len(cuts) else cuts[-1].mid() * 2 - 1
            mu = (lo + hi) / 2
            weight = {label: mu, other: big}
            f = lambda m: energy[m] + weight["+" if sector[m] == 1 else "-"] * gap[m]
            assert all(f(m) > f(v) for m in ns if m != v), (label, v)
            attainable.append(v)
    assert len(set(attainable)) == len(attainable) == 33
    print("attainable winners of the two-parameter family:", len(attainable))

    # The anti-diagonal mu+ = c, mu- = 1-c: affine combinations of the two
    # fixed-parity ground energies, (1-c) E_plus + c E_minus.
    anti = {n: (sector[n] * gap[n], energy[n] + (gap[n] if sector[n] == -1 else 0)) for n in ns}
    # F = E_plus + c * s(n) g(n); the chain runs from c=+inf to c=-inf.
    achain, abreaks = envelope(anti)
    assert achain == [10965, 11363, 11092, 11491, 11891, 11619, 11618, 11346, 11746, 11075,
                      11474, 11050, 11448, 11177, 11576, 11304, 11977, 11033, 11431]
    for c, winner in ((0, 11177), (1, 11618)):
        assert all(value(anti, m, c) > value(anti, winner, c) for m in ns if m != winner)
    print("fixed-parity interpolation winners, weight on E_minus from +inf to -inf:", achain)
    print("  breakpoints", [b.str(6) for b in abreaks])

    # --- monotone no-go ----------------------------------------------------------
    witnesses = {10875: (11177, 11681), 11737: (11618, 12003), 11488: (11618, 11346)}
    for n, (up, down) in witnesses.items():
        assert sector[up] == sector[down] == sector[n]
        assert energy[up] < energy[n] and gap[up] > gap[n]
        assert energy[down] < energy[n] and gap[down] < gap[n]
    print("10875, 11737, 11488 dominated in both gap directions inside their own sector")

    # Exact sector Pareto sets (vacuum energy down, gap up or down) from rigorous
    # float enclosures; every pair is decided.
    e_lo, e_hi = interval_bounds([energy[n] for n in ns])
    g_lo, g_hi = interval_bounds([gap[n] for n in ns])
    sig = np.array([sector[n] for n in ns])

    def front(idx, direction):
        """Members of idx not strictly dominated (lower energy, gap better in direction)."""
        kept = []
        for i in idx:
            lower_e = e_hi[idx] < e_lo[i]
            higher_e = e_lo[idx] > e_hi[i]
            if direction == 1:
                better_g, worse_g = g_lo[idx] > g_hi[i], g_hi[idx] < g_lo[i]
            else:
                better_g, worse_g = g_hi[idx] < g_lo[i], g_lo[idx] > g_hi[i]
            others = idx != i
            dominated = np.any(lower_e & better_g & others)
            # Not dominated only when every rival is strictly worse in one coordinate.
            assert dominated or np.all(higher_e | worse_g | ~others)
            if not dominated:
                kept.append(ns[i])
        return kept

    selectable, sizes = set(), {}
    for s in (1, -1):
        for direction in (1, -1):
            kept = front(np.nonzero(sig == s)[0], direction)
            sizes[(s, direction)] = len(kept)
            selectable |= set(kept)
    assert sizes == {(1, 1): 11, (1, -1): 40, (-1, 1): 10, (-1, -1): 33}
    assert len(selectable) == 92 and not selectable & {10875, 11737, 11488}
    assert {11431, 11618, 11177} <= selectable
    blind_up = front(np.arange(len(ns)), 1)
    blind_down = front(np.arange(len(ns)), -1)
    assert blind_up == [11092, 11363, 11431, 11491, 11618, 11619, 11763, 11891, 11892, 12036]
    blind = set(blind_up) | set(blind_down)
    assert len(blind_down) == 35 and len(blind) == 44 and blind <= selectable
    print("sector-monotone selectable zeros:", len(selectable), sizes)
    print("parity-blind monotone selectable zeros:", len(blind), "; gap-rewarding front:", blind_up)

    # Exact ranks of the five named zeros in the two native orders.
    ranks = {}
    for n in (10875, 11431, 11488, 11618, 11737):
        i = ns.index(n)
        below = e_hi < e_lo[i]
        assert np.all(below | (e_lo > e_hi[i]) | (np.arange(len(ns)) == i))
        above = g_lo > g_hi[i]
        assert np.all(above | (g_hi < g_lo[i]) | (np.arange(len(ns)) == i))
        ranks[n] = (1 + int(below.sum()), 1 + int(above.sum()))
    assert ranks == {10875: (1172, 588), 11431: (173, 1), 11488: (194, 223),
                     11618: (1, 158), 11737: (426, 1076)}
    print("(vacuum-energy rank, gap rank):", ranks)

    # --- minimal unification cost -----------------------------------------------
    e_min, g_max = energy[11618], gap[11431]
    assert all(energy[m] > e_min for m in ns if m != 11618)
    assert all(gap[m] < g_max for m in ns if m != 11431)

    def cost(n):
        a, b = energy[n] - e_min, g_max - gap[n]
        return (a + b + abs(a - b)) / 2

    costs = {n: cost(n) for n in ns}
    best = 11491
    assert all(costs[m] > costs[best] for m in ns if m != best)
    assert arb("0.0021434") < costs[best] < arb("0.0021439")
    assert costs[11763] - costs[best] > arb("0.0000335")
    l1 = energy[best] - e_min + g_max - gap[best]
    assert arb("0.0031699") < l1 < arb("0.0031700")
    assert arb("0.09055") < costs[best] / g_max < arb("0.09057")
    for n, (lo, hi) in {10875: ("0.026143", "0.026145"), 11488: ("0.0084160", "0.0084163"),
                        11737: ("0.0219155", "0.0219158")}.items():
        assert arb(lo) < costs[n] < arb(hi)
    print("minimal L-infinity unification cost at 11491:", costs[best].str(8),
          "runner-up 11763 margin >", "0.0000335")
    print("L1 cost at 11491:", l1.str(8), "  cost / max gap:", (costs[best] / g_max).str(6))


if __name__ == "__main__":
    verify()
