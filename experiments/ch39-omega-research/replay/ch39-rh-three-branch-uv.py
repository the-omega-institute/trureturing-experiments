"""Compare the three finite clock solutions above their respective crossings."""
from pathlib import Path
from runpy import run_path

from flint import arb, ctx

fixed = run_path(str(Path(__file__).with_name("ch39-rh-fixed-clock-crossings.py")))
positive = fixed["positive"]


def inverse(energy, mf, masses, d):
    e = sum(n*positive((energy/m).log()) for m, n in masses)
    t = (energy/1000).log()
    u = positive((energy/(mf*arb(8).sqrt())).log())
    v = positive((energy/(mf*arb(24).sqrt())).log())
    return (d[0]-(41*t+e)/10-(16*u+4*v)/5,
            d[1]+(19*t-e)/6, d[2]+7*t-2*u)


def window(location, d, left, right):
    q0, mf0, spec0, _ = left
    q1, mf1, spec1, _ = right
    common = (d[2]+7*(q0.union(q1)/1000).log()
              -2*positive(location-arb(8).log()/2))/(2*arb.pi())
    assert common > 1 and sum(n for m,n in spec0 if m < q1)>19
    print("Common inverse alpha:", common.str(8))
    outputs = []
    for cut in (2*arb.pi(), arb(0)):
        endpoints = []
        for mf, spec in [(mf1, spec1), (mf0, spec0)]:
            lo, hi = 10**14, 10**16
            assert lo > q0 and inverse(arb(lo), mf, spec, d)[1] > cut
            assert inverse(arb(hi), mf, spec, d)[1] < cut
            for _ in range(52):
                if hi-lo <= 1:
                    break
                mid = (lo+hi)//2
                sign = inverse(arb(mid), mf, spec, d)[1]-cut
                if sign > 0:
                    lo = mid
                elif sign < 0:
                    hi = mid
                else:
                    break
            endpoints.append(arb(lo).union(arb(hi)))
        band = endpoints[0].union(endpoints[1])
        assert inverse(band.upper(), mf1, spec1, d)[0] > cut
        outputs.append(band)
    print("First alpha = 1 and inverse zero (weak channel), GeV:", *(x.str(9) for x in outputs))
    return outputs


def main():
    ctx.prec, ctx.cap = 1536, 3
    data = fixed["inputs"]()
    d = [8*arb.pi()**2/arb(g)**2 for g in (".467453", ".63811", "1.0583")]
    ranges = []
    for lo, hi, sign in [(".31499",".31502",1), ("1.64484","1.64488",-1),
                         ("8.19470","8.19480",1)]:
        w0, w1 = arb(lo), arb(hi)
        left, right = [fixed["endpoint"](w, data) for w in (w0, w1)]
        assert sign*left[-1] < 0 and sign*right[-1] > 0
        ranges.append(window(w0.union(w1), d, left, right))
    for values in zip(*ranges):
        assert values[0] > values[1] and values[1] > values[2]


if __name__ == "__main__":
    main()
