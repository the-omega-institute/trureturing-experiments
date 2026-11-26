"""Both triangle orbits of the two-simplex cluster, with an exterior holonomy fixed."""
from fractions import Fraction as F
from functools import cache
from itertools import product
from pathlib import Path
from runpy import run_path

old = run_path(str(Path(__file__).with_name('ch39-weyl-glued-simplices.py')))
CHARGES, allowed, pack = old['CHARGES'], old['allowed'], old['pack']
# Use the same fifteen slots. Replace the old 014 coefficient by that of 012.
probe = (1, 4)
slot = old['FACES'].index((probe,))
FACES = old['FACES'][:slot]+(((1, 2),),)+old['FACES'][slot+1:]


@cache
def contract(d, types):
    return old['contract'](d, types, probe, FACES)


def coefficients():
    parts = [(*labels[:3], *pack(labels)) for labels in old['caps']()]
    cores = {vals: pack(vals) for vals in product(CHARGES, repeat=3)}
    rows = [[0]*(1 << 15) for _ in range(4)]
    # Edge 14 is measured. Edge 12, now integrated, fixes a monogon instead.
    for a, b, c in product(CHARGES, repeat=3):
        v, w = c-a, -b-c
        if v not in CHARGES or w not in CHARGES or not allowed((v, a, -c)) or not allowed((w, b, c)):
            continue
        left = {charge: pack((a, b, c, charge, v, w)) for charge in CHARGES}
        for d, e, f, rc, rw, rm in parts:
            for g in CHARGES:
                charge = -g-a-d
                if charge not in left or not allowed((charge, g, a, d)):
                    continue
                u, v0 = g-b-e, -g-c-f
                core = cores.get((u, v0, g))
                if core is None or not allowed((u, -g, b, e)) or not allowed((v0, g, c, f)):
                    continue
                cc, cw, cm = core
                lc, lw, lm = left[charge]
                mz, mc = contract(3, cc+lc+rc)
                wz, mw = contract(2, cw+lw+rw)
                z = mz*wz
                mask = cm+(lm << 3)+(rm << 9)
                for row, value in zip(rows, (z, z*(a+b)**2, 15*mz*mw, 15*wz*mc)):
                    row[mask] += value
    return rows


def geometric_face(face):
    vertex_set = {i for edge in face for i in edge}
    if len(face) == 1:
        vertex_set.add(0)
    return frozenset(vertex_set)


def verify(rows):
    import sympy as sp
    one = old['base']['face_coefficients']()
    faces = list(map(geometric_face, FACES))
    swap = {2: 4, 4: 2}
    indices = [faces.index({swap.get(i, i) for i in geometric_face(face)}) for face in old['base']['FACES']]
    units = [old['SCALE']**2]+[60*old['SCALE']**2]*3
    for subset in range(512):
        mask = sum(1 << j for i, j in enumerate(indices) if subset & (1 << i))
        assert all(F(row[mask], scale) == ref[subset] for row, scale, ref in zip(rows, units, one, strict=True))
    u = [y-w for y, w in zip(rows[1], rows[2])]
    v = [w-c for w, c in zip(rows[2], rows[3])]
    for coeff in (v, [5*x-112*y for x, y in zip(u, v)], [1006712*y-3025*x for x, y in zip(u, v)]):
        scores = [c*12**(15-i.bit_count()) for i, c in enumerate(coeff)]
        for j in range(15):
            bit = 1 << j
            for i in range(1 << 15):
                if i & bit:
                    scores[i] += scores[i ^ bit]
        assert min(scores) >= 0 and max(scores) > 0
    for active, ratio in [(0xc8, F(112, 5)), (0x676c, F(1006712, 3025))]:
        assert v[active] > 0 and u[active] == ratio*v[active]
        assert all(not (u[i] or v[i]) for i in range(len(v)) if i != active and i & active == i)
    poly = old['homogeneous_polynomials']
    positive = old['positive_on_interval']
    single = poly(one, 9, [1]*4)
    central = poly(old['coefficients'](), 15, units)
    z, *curves = poly(rows, 15, units)
    for sign, baseline, degree, leading in [
        (1, single, 5, (F(409, 162), F(7, 54), F(17, 324))),
        (-1, central, 3, (F(71, 45), F(1, 3), F(5, 18))),
    ]:
        base_z, *base_h = baseline
        for h, base_curve, first in zip(curves, base_h, leading, strict=True):
            diff = sign*(h*base_z-base_curve*z)
            assert all(diff.nth(i) == 0 for i in range(degree)) and diff.nth(degree) == first
            assert positive(diff)
    t = sp.Symbol('t')
    strengths = [4*z+sp.Poly(1+30*t,t)*p.exquo(sp.Poly(t,t)) for p in curves]
    a, _, c = strengths
    assert positive(a.diff()*c-c.diff()*a)
    lo, hi = sp.Rational('0.076393202250020'), sp.Rational('0.076393202250022')
    assert (3-10*hi)**2 < 5 < (3-10*lo)**2
    assert sp.Rational('1.01271545265499') < a.eval(lo)/c.eval(lo)
    assert a.eval(hi)/c.eval(hi) < sp.Rational('1.01271545265502')
    print('Exterior: same sharp slopes 112/5 and 1006712/3025. Equal-strength spread supremum 1.0127154526550.')
    print('All three stiffnesses increase from isolated to exterior to shared. First corrections at orders 5 and 3.')


if __name__ == '__main__':
    verify(coefficients())
