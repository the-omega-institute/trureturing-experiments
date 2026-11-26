"""How much weight the six high zeros need to cross the arithmetic pivot threshold."""
from pathlib import Path
from runpy import run_path

from flint import acb, arb, arb_mat, arb_series, ctx

root = Path(__file__).parent
inputs = run_path(str(root/"ch39-arithmetic-height-bandwidth.py"))
threshold = "266.6444299"


def gram(readings, ids):
    return arb_mat([[(i*readings[i]+j*readings[j])/(i+j) if i+j else arb(1)
                     for j in ids] for i in ids])


def pivot_certificate(outside, band, ids):
    o,b = gram(outside,ids),gram(band,ids)
    n = len(ids)
    def evaluate(y):
        mat = o+y*b
        leading = arb_mat([[mat[i,j] for j in range(n-1)] for i in range(n-1)])
        return mat[n-1,n-1]*leading.det()-arb(threshold)**2*mat.det()
    vandermonde = arb_mat([[arb(i)**j for j in range(n+1)] for i in range(n+1)])
    coefficients = vandermonde.solve(arb_mat([[evaluate(arb(i))] for i in range(n+1)]))
    assert all(coefficients[i,0]>0 for i in range(1,n+1))
    return evaluate, coefficients[0,0]


def capped_gain(count, cap):
    quotient, remainder = divmod(count, cap)
    return quotient*cap**2 + remainder**2


def verify():
    ctx.prec,ctx.cap = 3072,3
    times = inputs["times"]()
    s = arb_series([arb(1)/2,1],prec=3)
    xi = s*(s-1)/2*((-s/2)*arb.pi().log()).exp()*(s/2).gamma()*s.zeta()
    mass0 = xi[2]/xi[0]
    heights = [acb.zeta_zero(n).imag for n in range(11485,11491)]
    weight = sum(1/g**2 for g in heights)/mass0
    band = [sum(1/(g*g+j*j) for g in heights)/(weight*mass0) for j in range(len(times))]
    outside = [(t-weight*b)/(1-weight) for t,b in zip(times,band,strict=True)]
    for ids in ([0,1,2,3], [*range(8),28]):
        _, constant = pivot_certificate(outside,band,ids)
        assert constant>0
    evaluate, constant = pivot_certificate(outside,band,[0,1,2])
    assert constant<0
    lo,hi = arb(0),arb(1)
    for _ in range(120):
        mid = (lo+hi)/2
        value = evaluate(mid)
        assert value>0 or value<0
        if value>0:
            hi=mid
        else:
            lo=mid
    odds = lo.union(hi)
    fraction = odds/(1+odds)
    gain = odds*(1-weight)/weight
    assert arb("7982.06671") < gain < arb("7982.06672")
    assert arb(".01643732") < fraction < arb(".01643733")
    assert arb("2.09369e-6") < weight < arb("2.09370e-6")
    first_shift = (fraction-weight)*(band[1]-outside[1])
    assert arb("2.63397e-5") < first_shift < arb("2.63398e-5")
    assert 1/arb(2)**16 < first_shift < 1/arb(2)**15
    # Reconstruct the geometric fibres and rotation cycles using their own certificate.
    run_path(str(root/"ch39-intrinsic-shell-comparisons.py"))["verify"]()
    assert 160*3+16*9 == 624 and 160*3**2+16*9**2 == 2736 < gain
    assert capped_gain(624,12) < gain < capped_gain(624,13)
    for cap,expected in [(14,7872),(15,8544)]:
        value = 12*(capped_gain(16,cap)+capped_gain(36,cap))
        assert value == expected
        assert (value<gain) if cap==14 else (value>gain)
    assert 16*(9+3+3)+25*(5*3)+3*3 == 624
    assert 41*15**2+9**2 == 9306 > gain
    print("Original band weight",weight.str(19),"critical fraction",fraction.str(19),"gain",gain.str(19))
    print("Channel caps: arbitrary groups 13; whole cycles 15; separate fibres 15.")


if __name__ == "__main__":
    verify()
