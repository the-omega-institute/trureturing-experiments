"""Keep one pair of shell sites for the full record and identify its port law."""
from fractions import Fraction
from itertools import combinations
from pathlib import Path
from runpy import run_path

from flint import arb, fmpz_mat

model = run_path(str(Path(__file__).with_name("ch39-shell-selected-short-words.py")))
short, local = model["short"],model["local"]
events = short["events"]


def pole_bound():
    lap, pole = local["lifted_laplacian"]()
    inside = [i for i in range(25) if i != pole]
    whole = lap.charpoly()
    minor = fmpz_mat([[lap[i,j] for j in inside] for i in inside]).charpoly()
    common = whole.gcd(minor)
    num, den = minor//common, whole//common
    roots = whole.complex_roots()
    assert all(z.imag.is_zero() for z,_ in roots)
    eigenspaces = sorted(((z.real,m) for z,m in roots),key=lambda r:float(r[0].mid()))
    largest, mult = eigenspaces[-1]
    assert mult==1 and all(z<largest for z,m in eigenspaces[:-1])
    weight = num(largest)/den.derivative()(largest)
    assert weight > arb(".94476")
    pair_amplitude = 2*weight**2-1
    assert pair_amplitude>0 and pair_amplitude**2 > arb(".616")
    print("double-pole return probability at least", (pair_amplitude**2).lower().str(10))


def contraction_ranks(totals):
    rows = [tuple(x*x for x in word) for word in totals]
    assert all(fmpz_mat(four).det()!=0 for four in combinations(rows,4))
    # The cubic moment contracted with a normal to three points should have
    # column space spanned by the two remaining probability vectors.
    for triple in combinations(range(5),3):
        normal = [(-1)**j*fmpz_mat([[rows[i][k] for k in range(4) if k!=j]
                                    for i in triple]).det() for j in range(4)]
        dots = [sum(a*b for a,b in zip(normal,row)) for row in rows]
        assert [i for i,dot in enumerate(dots) if dot==0]==list(triple)
        matrix = fmpz_mat(4,4)
        for weight,dot,row in zip(totals.values(),dots,rows,strict=True):
            column = fmpz_mat([[k] for k in row])
            matrix += weight*dot*column*column.transpose()
        assert matrix.rank()==2
        complement = [rows[i] for i in range(5) if i not in triple]
        assert fmpz_mat(matrix.tolist()+complement).rank()==2


def verify():
    _,counts = model["radial_counts"]()
    words = model["selected"]
    totals = {v:sum(c for w,c in zip(words,counts) if w==v) for v in words}
    assert sorted(totals.values()) == [1,18,81,255,270]
    # All words have the same integer mass 5^17. Mix the complete record
    # probabilities, with the same denominator for every word at each length.
    for size in range(1,34):
        weights = {}
        for v,count in totals.items():
            current = events["orbit_weights"](size,v)
            for key,n in current.items():
                weights[key] = weights.get(key,0)+count*n
        assert sum(weights.values())==625*5**(17*size)
        if size<33:
            assert not events["possible_light_events"](weights)
    masses,_ = events["event_masses"](weights)
    masses["ref"] -= weights[18,6,3]
    ys = [Fraction(masses[n],masses["ref"]) for n in events["names"]]
    assert all(events["low"][n]<y<events["high"][n] for n,y in zip(events["names"],ys))
    ports = [sum(c*v[i]**2 for v,c in totals.items()) for i in range(4)]
    state = [arb(v)/ports[0] for v in ports[1:]]+[short["q"](y) for y in ys]
    errors = short["residuals"](state,short["q"]("26.81668495751949"))
    assert max(e.upper() for e in errors) < arb("3.06e-7")
    probs = [tuple(Fraction(x*x,5**17) for x in w) for w in totals]
    contraction_ranks(totals)
    # Identical-port coincidences distinguish this mixture from fresh draws.
    mean = [Fraction(w,625*5**17) for w in ports]
    excess = sum(Fraction(c,625)*sum((a-b)**2 for a,b in zip(p,mean)) for c,p in zip(totals.values(),probs))
    assert excess>0
    pole_bound()
    print("fixed shell pair 33-copy gauge residuals",*(e.upper().str(9) for e in errors))
    print("excess equal-port coincidence",short["q"](excess).str(9),"five points in general position")


if __name__ == "__main__":
    verify()
