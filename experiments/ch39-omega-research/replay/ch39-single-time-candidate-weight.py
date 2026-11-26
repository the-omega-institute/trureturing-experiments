"""Optimal one-time upper bounds on the six-zero band's weight at unit scale."""
from pathlib import Path
from runpy import run_path

from flint import acb, arb, arb_series, ctx

root = Path(__file__).parent
pivot = run_path(str(root/"ch39-candidate-band-pivot-gain.py"))


def verify():
    ctx.prec,ctx.cap = 240,3
    half = arb(1)/2
    s = arb_series([half,1],prec=3)
    xi = s*(s-1)/2*((-s/2)*arb.pi().log()).exp()*(s/2).gamma()*s.zeta()
    mass0 = xi[2]/xi[0]
    first = acb.zeta_zero(1).imag**2
    squares = [acb.zeta_zero(n).imag**2 for n in range(11485,11491)]
    mass = sum(1/g for g in squares)
    logpi = arb.pi().log()
    def reading(j):
        s = j+half
        if j<120:
            z = arb_series([s,1],prec=2).zeta()
            derivative = z[1]/z[0]
        else:
            # Sum log(n)/n^s for n>=2, bounded by its first term plus an integral.
            derivative = arb(0, arb(2)**(1-j))
        return (1/s+1/(s-1)+((s/2).digamma()-logpi)/2+derivative)/(2*j*mass0)
    def band(j):
        return sum(1/(g+j*j) for g in squares)/mass
    t = [arb(1),reading(1),reading(2)]
    w = mass/mass0
    b = [band(i) for i in range(3)]
    o = [(a-w*y)/(1-w) for a,y in zip(t,b)]
    threshold, _ = pivot["pivot_certificate"](o,b,[0,1,2])
    f0,f1 = arb(".01643732111729698"),arb(".01643732111729700")
    assert threshold(f0/(1-f0)) < 0 < threshold(f1/(1-f1))
    critical = f0.union(f1)
    def values(j):
        t = reading(j)
        base = first/(first+j*j)
        bands = (band(j),squares[0]/(squares[0]+j*j))
        assert base < t < bands[1]
        return [((t-base)/(v-base),base+critical*(v-base)-t) for v in bands]
    min_idx = (9708,9707)
    min_values = [values(j)[k][0] for k,j in enumerate(min_idx)]
    best = values(8709)
    windows = [(5600,16386),(5601,16378)]
    cutoff = 40000
    for j in range(1,cutoff+1):
        for k,(upper, margin) in enumerate(values(j)):
            if j != min_idx[k]:
                assert upper > min_values[k]
            if j != 8709:
                assert margin < best[k][1]
            lo,hi = windows[k]
            assert margin>0 if lo<=j<=hi else margin<0
    tail = reading(cutoff)
    assert (cutoff**2*tail-first)/max(squares) > critical
    assert 1-tail > critical
    assert arb(".01440186") < min_values[1] < arb(".01440188")
    assert arb(".00121278") < best[1][1] < arb(".00121279")
    print("Minimum upper bounds",[(j,v.str(16)) for j,v in zip(min_idx,min_values)])
    print("Largest absolute readout gaps at 8709",[g.str(16) for _,g in best])
    print("Complete excluding offset intervals",windows)


if __name__ == "__main__":
    verify()
