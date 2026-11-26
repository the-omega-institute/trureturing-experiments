"""First negative color harmonic of the Newton power and a nearby integer microscopic matching."""
from fractions import Fraction as F
from pathlib import Path
from runpy import run_path
from flint import acb, arb

ROOT = Path(__file__).resolve().parent
high = run_path(str(ROOT/"ch39-high-zero-planck-window.py"))
model, q = high["model"], high["q"]
L0 = high["L"]
N = 30484476
K = N+2


def match_at_integers(prefactor):
    steps = int(64*L0)
    bounds = [(q(c*F(999,1000)),q(c*F(1001,1000))) for c in model["GAUGE"]]
    full = [a.union(b) for a,b in bounds]
    def at_shift(gauge):
        derivatives = []
        def record(j, inv, deriv):
            if j in (steps-1,steps):
                derivatives.extend(deriv)
        _, endpoint = model["running"]["integrate"](gauge,steps)
        model["running"]["integrate"](gauge,steps+1,record)
        assert derivatives and all(abs(x)<1 for x in derivatives)
        return [v+arb(0,q("0.000000016")) for v in endpoint]
    ratios = [q(v) for v in high["V"]]
    full_values = at_shift(full)
    outer_low = max((x/y).lower() for x,y in zip(full_values,ratios))
    outer_high = min((x/y).upper() for x,y in zip(full_values,ratios))
    low, upper = [], []
    for i in range(3):
        for side, dest in [(0,upper),(1,low)]:
            face = full.copy()
            face[i] = bounds[i][side]
            x = at_shift(face)[i]/ratios[i]
            dest.append(x.lower() if side==0 else x.upper())
    for power in (N,N+1):
        shift = (prefactor/power).log()/2
        assert abs(shift) < q("1.6e-8")
        coeff = 4*power*(3*arb(5).sqrt()/2)
        def labels(left,right):
            counts = [(coeff/v-q(high["OFFSET"])).sqrt().zeta_nzeros().unique_fmpz() for v in (right,left)]
            assert all(x is not None for x in counts)
            return list(range(int(counts[0])+1,int(counts[1])+1))
        assert labels(outer_low,outer_high) == labels(max(low),min(upper)) == list(range(11485,11491))
        print("Integer action power and exact-Newton log-scale shift:",power,shift.str(12))


def verify():
    prefactor = (high["geo"]["BOOK_PLANCK"]/(1000*q(L0).exp()))**2/(8*arb.pi())
    theta = prefactor-N
    assert q(".09") < theta < q(".11")
    assert theta*(1-theta)/(K*(K-1)) > q(F(1,2**54))
    for j in range(11485,11491):
        height = acb.zeta_zero(j).imag
        amplitude = 16/(height**2+q(high["OFFSET"]))
        assert q(F(1,2**24)) < amplitude/2 < q(F(1,2**23))
        assert (1-2*amplitude)**(-prefactor-2) < 2**13
    e = F(1,2**26)
    assert K*K*e*e < 2
    x = q("2.56e-8")
    point = acb((1-x*x)/(1+x*x),2*x/(1+x*x))
    assert (point**K).real > 0 and 2*K*x < 2 and abs(point-1) > q("5.1e-8")
    match_at_integers(prefactor)


if __name__ == "__main__":
    verify()
