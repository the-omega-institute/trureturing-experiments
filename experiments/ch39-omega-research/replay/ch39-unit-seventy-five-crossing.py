"""一份七十五维重物质在同一十九点连接及三辐边修正下的平滑三线交点。"""
from collections import defaultdict
from fractions import Fraction
from itertools import combinations
from pathlib import Path
from runpy import run_path

from flint import acb, arb, ctx, fmpz_mat

smooth = run_path(str(Path(__file__).with_name("ch39-smooth-two-stiffness-crossing.py")))
matter, soft = smooth["matter"], smooth["soft"]
beta = smooth["beta"]


def index_table():
    # End(wedge^2 5) = 1+24+75, End(5) = 1+24.
    rows = defaultdict(lambda: [Fraction(0) for _ in range(3)])
    for degree, sign in ((2, 1), (1, -1)):
        weights = list(combinations(range(5), degree))
        for left in weights:
            for right in weights:
                v = [int(i in left)-int(i in right) for i in range(5)]
                q = -2*sum(v[:3])+3*sum(v[3:])
                terms = (Fraction(q*q, 60), Fraction((v[3]-v[4])**2, 4),
                         Fraction((v[0]-v[1])**2, 4))
                for j,t in enumerate(terms):
                    rows[abs(q)][j] += sign*t
    assert dict(rows) == {0:[0,16,12], 5:[15,9,12], 10:[10,0,1]}
    return {q:[int(t) for t in row] for q,row in rows.items()}


def shell_spectra(charges=(0, 5, 10), zero_index=11431):
    height = acb.zeta_zero(zero_index).imag
    base, neutral, _ = soft["vector_pencil"](height)
    matrix = base+neutral
    integer = fmpz_mat([[int(t.real.unique_fmpz()) for t in row] for row in matrix.tolist()])
    neutral_mass = [z.real for z,k in integer.charpoly().complex_roots() if z!=0 for _ in range(k)]
    # All three radial edges have unit transport for every charge. The norm
    # of the rank-three stiffness removal equals four times its per-edge weight.
    amount = arb(23)/25*arb(2)**-20/(height**2+soft["point"]["study"]["offset"])
    shift = arb(0).union(4*amount)
    values = {q: smooth["masses"](base+matter["nonspokes"](q, height))
              for q in charges if q}
    if 0 in charges:
        assert len(neutral_mass) == 18
        values[0] = neutral_mass
    result = {q:[v-shift for v in masses] for q,masses in values.items()}
    assert all(v>arb(".01") for row in result.values() for v in row)
    return result


def verify_decrease(masses, indices, radius, box):
    # d sheet / d log(momentum) integrates z/(m2+z) with positive weight,
    # z=(radius*f)^2*(1-y*y)/4. Positive (vector) terms need upper bounds.
    def slope(levels, factor, vector):
        if factor is None:
            return len(levels)
        def primitive(y):
            return (8*y-y**3)/7 if vector else (3*y-y**3)/2
        value = arb(0)
        for i in range(32):
            lo, hi = arb(i)/32, arb(i+1)/32
            y = lo if vector else hi
            z = (radius*factor)**2*(1-y*y)/4
            value += (primitive(hi)-primitive(lo))*sum(z/(v+z) for v in levels)
        return value

    matter_levels = {q: [v+box for v in (row if q else [arb(0), *row])]
                     for q, row in masses.items()}
    intervals = [(arb(i)/4, arb(i+1)/4) for i in range(4,8)] + [(arb(2), None)]
    bounds = []
    for lo, hi in intervals:
        neutral, heavy = (slope(masses[q], hi, True) for q in (0, 5))
        dirac = {q: slope(row, lo, False) for q, row in matter_levels.items()}
        for i, (b, n, h) in enumerate(zip(beta, (0,2,3), (5,3,2))):
            upper = (-b/2+7*(n*neutral+h*heavy)/4
                     -2*sum(indices[q][i]*v for q,v in dirac.items())/3)/arb.pi()
            assert upper < -arb("21.37")
            bounds.append(upper.upper())
    print("Upper bound on all log-momentum slopes above matching", max(bounds).str(9))


def verify():
    ctx.prec = 224
    indices = index_table()
    mass = shell_spectra()
    geo = soft["geo"]
    radius = geo["theta"](804)/geo["theta"](12)
    mz, em, weak = arb("91.1876"), arb("127.930"), arb(".23129")
    initial = [3*(1-weak)*em/5, em*weak, 1/arb(".1180")]
    pi = arb.pi()
    def intercept(x, factor=1):
        q2 = (radius*factor)**2
        vectors = [sum(smooth["sheet"](v,q2,True) for v in mass[q]) for q in (0,5)]
        def dirac(q, m2):
            levels = mass[q] if q else [arb(0), *mass[0]]
            return sum(smooth["sheet"](v+m2,q2,False) for v in levels)
        exposures = {q:dirac(q,x.upper()).union(dirac(q,x.lower())) for q in indices}
        return [a+7*(n*vectors[0]+b*vectors[1])/(4*pi)
                -2*sum(indices[q][i]*v for q,v in exposures.items())/(3*pi)
                for i,(a,n,b) in enumerate(zip(initial,(0,2,3),(5,3,2)))]
    def differences(values):
        return values[0]-values[2], values[1]-values[2]
    def projection(values):
        first, second = differences(values)
        return (beta[1]-beta[2])*first-(beta[0]-beta[2])*second

    # Locate the mass from the component independent of the light running.
    lo, hi = arb(".01"), arb(".1")
    assert projection(intercept(lo))*projection(intercept(hi)) < 0
    positive_left = projection(intercept(lo)) > 0
    for _ in range(36):
        mid = (lo+hi)/2
        sign = projection(intercept(mid))
        assert sign > 0 or sign < 0
        if (sign > 0) == positive_left:
            lo = mid
        else:
            hi = mid
    box = lo.union(hi)+arb(0,arb("1e-9"))
    L = 2*pi*differences(intercept(lo.union(hi)))[0]/(beta[0]-beta[2])
    logs = L+arb(0,arb("1e-5"))
    assert box > arb(".01")
    reference = radius*(-logs).exp()
    error = 19*reference**2/(8*arb(".01"))*(35/(4*pi)+50/(3*pi))
    def readings(x, scale, factor=1):
        return [a-b*(scale+arb(factor).log())/(2*pi)+arb(0,error.upper())
                for a,b in zip(intercept(x,factor),beta)]
    assert projection(readings(box.lower(),logs.mid()))*projection(readings(box.upper(),logs.mid())) < 0
    assert differences(readings(box,logs.lower()))[0] > 0
    assert differences(readings(box,logs.upper()))[0] < 0
    common = readings(box,logs)[2]
    assert arb("1.3680e16") < mz*logs.exp() < arb("1.3682e16")
    assert arb(".02251147") < box < arb(".02251149")
    assert arb("17.2031") < common < arb("17.2033")
    print("One Dirac 75, common shell, smooth matching scale GeV, mass squared, inverse",
          (mz*logs.exp()).str(13),box.str(12),common.str(12))
    later = readings(box, logs, 2)
    assert all(v < -arb("6.85") for v in later)
    print("At twice the matching momentum", *[v.str(10) for v in later])

    verify_decrease(mass, indices, radius, box)
    # Disjoint brackets certify the first zero above matching for every root
    # in the crossing rectangle. The subtraction momentum and masses stay fixed.
    for i, lower, upper in ((1, "1.679471", "1.679481"),
                            (0, "1.679495", "1.679504"),
                            (2, "1.682196", "1.682207")):
        assert readings(box, logs, arb(lower))[i] > 0
        assert readings(box, logs, arb(upper))[i] < 0
        print(("Y","W","C")[i], "zero momentum / matching momentum", lower, upper)


if __name__ == "__main__":
    verify()
