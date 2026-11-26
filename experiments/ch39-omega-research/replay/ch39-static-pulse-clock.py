"""Rank-one arithmetic mixing at one boundary of a seven-step propagation clock."""
from itertools import islice
from pathlib import Path
from runpy import run_path

from flint import acb, arb, fmpq, fmpz_poly

pulses = run_path(str(Path(__file__).with_name("ch39-finite-port-pulses.py")))


def free_ground(score, length=7):
    cosine = score.sqrt()
    if length != 7:
        order = length+1
        x = fmpz_poly([0,1])
        older, polynomial = fmpz_poly([1]), x
        for _ in range(1,order):
            older, polynomial = polynomial, 2*x*polynomial-older
        left = 1-fmpq(2,order**2)
        brackets = []
        for endpoint in (cosine.lower(),cosine.upper()):
            target = endpoint.fmpq()
            lo,hi = left,fmpq(1)
            assert polynomial(lo) <= target <= polynomial(hi)
            for _ in range(70):
                mid=(lo+hi)/2
                if polynomial(mid)<target:
                    lo=mid
                else:
                    hi=mid
            brackets.append(arb(lo).union(arb(hi)))
        return 2*(1-brackets[0].union(brackets[1]))
    for _ in range(3):
        cosine = ((1+cosine)/2).sqrt()
    return 2*(1-cosine)


def endpoint(energy, diagonal, length=7):
    pivot, derivative = diagonal-energy, arb(-1)
    for _ in range(length-1):
        derivative = -1+derivative/pivot**2
        pivot = 2-energy-1/pivot
    inverse = 1/(3-energy-1/pivot)
    return inverse, inverse**2*(1-derivative/pivot**2)


def verify():
    data = []
    for n, zero in enumerate(acb.zeta_zeros(10815, 1278), 10815):
        seq = pulses["pulse_records"](pulses["generator"](zero.imag))
        p = next(islice(seq, 6, 7))[0]
        data.append((n, 1/(arb(1)/4+zero.imag**2), p, free_ground(p)))
    total = sum(row[1] for row in data)
    data = [(n, w/total, p, e) for n, w, p, e in data]
    by_energy = sorted(data, key=lambda row: float(row[3].mid()))
    peak, runner = by_energy[:2]
    assert peak[0] == 11431 and runner[0] == 11831
    assert peak[3] < runner[3] and all(row[3] > runner[3] for row in by_energy[2:])

    def resolvent(energy, need_norm=False):
        a, da = endpoint(energy, 1)
        b, db = endpoint(energy, 3)
        value, norm, peak_norm = arb(0), arb(0), arb(0)
        for n, weight, score, _ in data:
            m = b+score*(a-b)
            denom = 1-2*m
            value += weight*m/denom
            if need_norm:
                term = weight*(db+score*(da-db))/denom**2
                norm += term
                if n == 11431:
                    peak_norm = term
        return value, norm, peak_norm

    def isolate(inverse_strength):
        lead, second = peak[3], runner[3]
        lo = (lead+arb("1e-7")*(second-lead)).mid()
        hi = (second-arb("1e-7")*(second-lead)).mid()
        assert resolvent(lo)[0]+inverse_strength < 0 < resolvent(hi)[0]+inverse_strength
        while hi-lo > arb("5e-12"):
            mid = ((lo+hi)/2).mid()
            test = resolvent(mid)[0]+inverse_strength
            if test < 0:
                lo = mid
            else:
                assert test > 0
                hi = mid
        return lo.union(hi)

    unit = isolate(1)
    top = isolate(0)
    _, norm, same_label = resolvent(unit, True)
    print("unit reference coupling: ground",unit.str(10),"label 11431",(same_label/norm).str(8))
    _, max_norm, _ = resolvent(top, True)
    _, da = endpoint(peak[3], 1)
    _, db = endpoint(peak[3], 3)
    residue = peak[1]/(4*(db+peak[2]*(da-db)))
    fraction = residue/((top-peak[3])**2 * max_norm)
    gap = runner[3]-top
    assert fraction > arb(".9936") and gap > arb(".00001937")
    assert arb(".9936") < same_label/norm < arb(".9938")
    print("all repulsive couplings: bare ground overlap >=", fraction.lower().str(9),
          "many-label gap >=",gap.lower().str(9))
    beta = 2**19
    # All other eigenvalues of any single 8-site clock exceed the half-angle edge.
    higher_floor = free_ground(arb(0))
    excited = sum((beta*(top-e)).exp() for _, _, _, e in by_energy[1:])
    excited += 1278*(8*52-1)*(beta*(top-higher_floor)).exp()
    probability = fraction/(1+excited)
    assert probability > arb(".9935")
    print("Gibbs label probability at beta >= 524288 exceeds", probability.lower().str(9))


if __name__ == "__main__":
    verify()
