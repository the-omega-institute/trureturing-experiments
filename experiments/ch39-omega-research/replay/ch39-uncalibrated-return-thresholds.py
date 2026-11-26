"""Identifiability of a zeta phase from old-site returns with unknown sector gains."""
from dataclasses import dataclass
from fractions import Fraction as F
from itertools import combinations
from pathlib import Path
from runpy import run_path

from flint import acb, arb, ctx

root = Path(__file__).parent
port = run_path(str(root/"ch39-cayley-port-readout.py"))
grid = run_path(str(root/"ch39-noisy-threshold-grid.py"))
shell, ball = port["shell"], port["ball"]


@dataclass(frozen=True)
class Profile:
    old_square: F
    bridge_square: F
    old_neutral: F
    old_charged: F
    mixed_neutral: F
    mixed_charged: F
    relay: F


def profiles(nodes, edges):
    neighbors = [{} for _ in nodes]
    for i, j, weight, charge in edges:
        neighbors[i][j] = weight, charge
        neighbors[j][i] = weight, -charge
    result = []
    for i, row in enumerate(neighbors[:19]):
        squares = [F(0), F(0)]
        triangles = [[F(0), F(0)] for _ in range(3)]
        for j, (weight, _) in row.items():
            squares[int(j >= 19)] += weight**2
        for j, k in combinations(row, 2):
            if k not in neighbors[j]:
                continue
            a, h = row[j]
            b, l = neighbors[j][k]
            c, m = neighbors[k][i]
            winding = abs(h+l+m)
            assert winding <= 1
            triangles[int(j >= 19)+int(k >= 19)][winding] += 2*a*b*c
        assert triangles[2][1] == 0
        result.append(Profile(*squares, *triangles[0], *triangles[1], triangles[2][0]))
    return result


def triple_difference(rows, attr):
    # Divide each third moment by its bridge square. The two nuisance columns
    # become 1 and relay/bridge_square, so an affine second difference kills both.
    r = [p.relay/p.bridge_square for p in rows]
    values = [getattr(p, attr)/p.bridge_square for p in rows]
    return (r[1]-r[2])*values[0]+(r[2]-r[0])*values[1]+(r[0]-r[1])*values[2]


def classify(nodes, data):
    groups = {}
    for node, profile in zip(nodes[:19], data, strict=True):
        key = tuple(sorted(node))
        assert groups.setdefault(key, profile) == profile
    assert len(groups) == 5
    pairs = list(combinations(groups.values(), 2))
    assert all(a.old_square*b.bridge_square != b.old_square*a.bridge_square for a, b in pairs)
    assert all(a.relay*b.bridge_square != b.relay*a.bridge_square for a, b in pairs)
    # At unit gains, a changed phase can be hidden by two changes proportional
    # to delta x. Bound those coefficients for every old-site pair.
    for a, b in pairs:
        def row(p):
            return p.relay/p.bridge_square, (p.old_charged+p.mixed_charged)/p.bridge_square
        r1, f1 = row(a)
        r2, f2 = row(b)
        drift = (f1-f2)/(r1-r2)
        offset = f1-r1*drift
        assert abs(drift) < 1 and abs(offset) < 1
        assert all(p.relay*drift+p.bridge_square*offset == p.old_charged+p.mixed_charged
                   for p in (a, b))
    exceptional = []
    for keys in combinations(sorted(groups), 3):
        rows = [groups[k] for k in keys]
        old = triple_difference(rows, "old_charged")
        mixed = triple_difference(rows, "mixed_charged")
        assert old or mixed
        if old and -mixed/old > 0:
            exceptional.append((keys, -mixed/old))
    assert [(tuple(p[0] for p in keys), int(1000*q)) for keys, q in exceptional] == [
        ((3, 21, 39), 427), ((21, 33, 39), 103695), ((33, 39, 51), 69)]
    print("Exceptional triples and positive squared gain ratios",
          [(k, float(q)) for k, q in exceptional])
    return exceptional


def recover(rows, seconds, thirds):
    pole, outer, active = rows
    bridge2 = seconds[0]/ball(pole.bridge_square)
    old2 = (seconds[1]-bridge2*ball(outer.bridge_square))/ball(outer.old_square)
    old = old2.sqrt()
    residuals, slopes = [], []
    for p, reading in zip(rows, thirds, strict=True):
        neutral = old*(old2*ball(p.old_neutral) + bridge2*ball(p.mixed_neutral))
        residuals.append((reading-neutral)/(bridge2*ball(p.bridge_square)))
        slopes.append(old*(old2*ball(p.old_charged) + bridge2*ball(p.mixed_charged)))
    r = [ball(p.relay/p.bridge_square) for p in rows]
    relay = (residuals[1]-residuals[0])/(r[1]-r[0])
    offset = residuals[0]-r[0]*relay
    assert slopes[0] == 0 and slopes[1] == 0
    x = (residuals[2]-r[2]*relay-offset)*bridge2*ball(active.bridge_square)/slopes[2]
    return x, (old, bridge2.sqrt(), relay, offset)


def observe(nodes, edges, indices, phase, gains, onsite):
    old, bridge, relay = gains
    scaled = [(i, j, (old if j < 19 else bridge if i < 19 else relay)*weight, charge)
              for i, j, weight, charge in edges]
    matrix = shell["matrices"](nodes, scaled, phase)
    for i in range(len(nodes)):
        matrix[i,i] = ball(onsite[int(i >= 19)])
    second = matrix**2
    third = second*matrix
    u = matrix[indices[0],indices[0]].real
    raw2 = [second[i,i].real for i in indices]
    raw3 = [third[i,i].real for i in indices]
    return u, raw2, raw3


def centered(u, raw2, raw3):
    return [v-u*u for v in raw2], [v-3*u*w+2*u**3 for w, v in zip(raw2, raw3, strict=True)]


def gain_from_threshold(height2, p):
    return (height2/((2*p/(1-p)**2)**2-arb(1)/4)).sqrt()


def threshold_pairings(heights, ratios):
    assert len(heights) == len(ratios) == 6
    ppm = ball(F(1, 10**6))
    stretch = (1+ppm)/(1-ppm)
    # Two compatible fits may be on opposite sides of a noisy measurement.
    uncertainty = (stretch**2).union(1/stretch**2)
    for i, h2 in enumerate(heights):
        for j, p in enumerate(ratios):
            gain = gain_from_threshold(h2, p*uncertainty)
            assert gain > 623 and gain < 625
            if i == j:
                assert gain.contains(624)
            else:
                assert abs(gain-624) > arb(".010")
    # With exact local heights, the earliest overlapping ratio bands have a
    # geometric-mean contact. Up to 2 ppm exclude all gains except 624.
    stretch = (1+2*ppm)/(1-2*ppm)
    for h2 in heights:
        for p in ratios:
            gain = gain_from_threshold(h2, p*(stretch**2).union(1/stretch**2))
            assert gain > 623 and gain < 625
    contacts = [(b.sqrt()-a.sqrt())/(b.sqrt()+a.sqrt()) for a,b in zip(ratios, ratios[1:])]
    assert all(t > contacts[-1] > 0 for t in contacts[:-1])
    assert arb("1.71272") < contacts[-1]/ppm < arb("1.71273")
    low, high = ratios[-2:]
    radius = contacts[-1]
    middle = (low*high).sqrt()
    assert (low*(1+radius)/(1-radius)).overlaps(middle)
    assert (high*(1-radius)/(1+radius)).overlaps(middle)
    print("36 pairings; first ambiguous ratio at ppm", (radius/ppm).str(15))


def verify():
    ctx.prec = 224
    nodes, edges = shell["geometry"]()
    data = profiles(nodes, edges)
    classify(nodes, data)
    indices = [nodes.index(p) for p in ((51,51,51), (3,15,87), (33,33,75))]
    rows = [data[i] for i in indices]
    assert rows[0].old_square == 0 and rows[2].old_charged > 0 and rows[2].mixed_charged > 0
    gains = (F(4,5), F(6,5), F(9,10))
    onsite = (F(1,100), F(-1,50))
    error = arb(0, "1e-24")
    factor = 1 + arb(0, "1e-6")
    cuts = [ball(F(q*q-1, q*q+1)) for q in port["cut_denominators"]]
    heights, ratios = [], []
    for n, model in grid["models"]():
        rho = acb.zeta_zero(n)
        phase = (rho-1)/rho
        u, raw2, raw3 = observe(nodes, edges, indices, phase, gains, onsite)
        x, recovered = recover(rows, *centered(u, raw2, raw3))
        for v, target in zip(recovered, (*gains, onsite[1]-onsite[0]), strict=True):
            assert v.overlaps(ball(target)) and v.rad() < arb("1e-45")
        assert x.overlaps(phase.real) and x.rad() < arb("1e-44")
        noisy = centered(u+error, [v+error for v in raw2], [v+error for v in raw3])
        x, _ = recover(rows, *noisy)
        assert all(x < cut or x > cut for cut in cuts)
        assert 11485 + sum(x > cut for cut in cuts) == n
        height2 = (1+x)/(4*(1-x))
        # Use independently formed E and D threshold positions from the existing
        # forward model, including separate measurement errors on the two locations.
        fermions, _, _ = grid["mass_ranges"](model, (F(0), F(0)), arb(1))
        p = fermions[0][0]/fermions[1][0]
        gain = gain_from_threshold(height2, p*factor/factor)
        assert gain > arb("623.98") and gain < arb("624.02") and gain.contains(624)
        assert gain_from_threshold(rho.imag**2, p).overlaps(arb(624))
        double = centered(u+2*error, [v+2*error for v in raw2], [v+2*error for v in raw3])
        other_x, _ = recover(rows, *double)
        heights.append((1+other_x)/(4*(1-other_x)))
        ratios.append(p)
        print(n, "height from noisy moments", height2.sqrt().str(12), "gain", gain.str(9))
    threshold_pairings(heights, ratios)


if __name__ == "__main__":
    verify()
