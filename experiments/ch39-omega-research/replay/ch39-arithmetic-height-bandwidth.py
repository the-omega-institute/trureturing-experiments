"""Sharp continuous height span forced by finitely many normalized RH time readings."""
import json
from fractions import Fraction as F
from math import prod
from pathlib import Path

from flint import acb, acb_mat, arb, arb_mat, arb_series, ctx

root = Path(__file__).parent


def ball(x):
    x = F(x)
    return arb(x.numerator)/x.denominator


def times():
    archived = json.loads((root/"ch39-minimum-tick-finite-inputs.json").read_text())["levels"]
    values = []
    for j in range(42):
        s = arb_series([arb(1)/2+j,1],prec=3)
        xi = s*(s-1)/2*((-s/2)*arb.pi().log()).exp()*(s/2).gamma()*s.zeta()
        v = xi[2]/xi[0] if j==0 else xi[1]/(2*j*xi[0])
        values.append(v)
    result = [v/values[0] for v in values]
    assert all((v*2**288).floor().unique_fmpz()==archived[j] for j,v in enumerate(result[1:],1))
    return result


def cleared_moments(readings, k):
    rates = [i*i for i in range(1,k+1)]
    coefficients = [1/prod(1-F(q,p) for q in rates if q!=p) for p in rates]
    moments = []
    for j in range(k+1):
        constant = ball(F(1,prod(rates))) if j==k else arb(0)
        value = constant + sum(ball(a*F(-1,p)**j)*v for a,p,v in
                               zip(coefficients,rates,readings[1:k+1],strict=True))
        moments.append(200**j*value)
    return moments


def gaussian_span(readings,k):
    m = cleared_moments(readings,k)
    size = (k+1)//2
    h = arb_mat([[m[i+j] for j in range(size)] for i in range(size)])
    g = arb_mat([[m[i+j+1] for j in range(size)] for i in range(size)])
    eig = acb_mat(h.solve(g)).eig()
    nodes = sorted((v.real for v in eig),key=lambda v:float(v.mid()))
    assert nodes[0]>0 and all(a<b for a,b in zip(nodes,nodes[1:]))
    vandermonde = arb_mat([[v**j for v in nodes] for j in range(size)])
    weights = vandermonde.solve(arb_mat([[v] for v in m[:size]]))
    assert all(weights[i,0]>0 for i in range(size))
    assert all(sum(weights[i,0]*v**j for i,v in enumerate(nodes)).overlaps(t) for j,t in enumerate(m))
    return (nodes[-1]/nodes[0]).sqrt()


def verify():
    ctx.prec,ctx.cap = 1536,3
    t = times()
    m0,m1,m2 = cleared_moments(t,2)
    k = m0*m2/m1**2
    ratio = k.sqrt()+(k-1).sqrt()
    assert arb("2.75415260") < ratio < arb("2.75415261")
    low_zero, high_zero = [acb.zeta_zero(n).imag for n in (11485,11490)]
    assert high_zero/low_zero < arb("1.00034")
    for j, bound in [(3520,high_zero/ratio),(36574,low_zero*ratio)]:
        assert acb.zeta_zero(j).imag < bound < acb.zeta_zero(j+1).imag
    lifts = [(3520,11490), (11485,36575)]
    for lo,hi in lifts:
        rr = (acb.zeta_zero(hi).imag/acb.zeta_zero(lo).imag)**2
        mean = m1/m0
        disc = ((1+rr)*mean)**2-4*rr*m2/m0
        for sg in (-1,1):
            a = ((1+rr)*mean+sg*disc.sqrt())/(2*rr)
            b = rr*a
            assert 0<a<mean<b
            q = (b-mean)/(b-a)
            w,v = [weight*m0*(1+node/200)*(1+4*node/200)
                   for weight,node in ((q,a),(1-q,b))]
            assert w>0 and v>0 and (w+v).contains(1)
            for j in (1,2):
                assert (w/(1+j*j*a/200)+v/(1+j*j*b/200)).overlaps(t[j])
    print("Two readings: minimum continuous height span",ratio.str(18), "boundary pairs",lifts)
    for count, lower in [(3,"3.34590"),(5,"5.57922"),(7,"7.95967"),(15,"20.12682"),(41,"86.32188")]:
        span = gaussian_span(t,count)
        assert arb(lower) < span < arb(lower)+arb(".00002")
        print(count,"readings: minimum continuous height span",span.str(18))


if __name__ == "__main__":
    verify()
