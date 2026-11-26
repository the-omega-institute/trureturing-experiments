"""由贝伦德壳的零点相位识别三十三阶局部展开的剩余系数。"""
from fractions import Fraction
from math import prod
from pathlib import Path
from runpy import run_path

from flint import acb, arb, ctx, fmpq_poly


def number(x):
    return arb(x.numerator)/x.denominator


def verify():
    geometry = run_path(str(Path(__file__).with_name("ch39-arc-relative-determinant.py")))
    ctx.prec = 224
    n = 18
    charges = [6*j for j in range(1,n+1)]
    v = [Fraction(1,prod(x*x-y*y for y in charges if y!=x)) for x in charges]
    step = min(Fraction(1,8*n*q*q)/abs(w) for q,w in zip(charges,v))/2
    coefficient = 2**n*step*v[-1]
    x = fmpq_poly([0,1])
    previous, chebyshev = fmpq_poly([1]), x
    direct = fmpq_poly()
    for w in v:
        term = 2*step*w
        direct += chebyshev*term.numerator/term.denominator
        previous, chebyshev = chebyshev, 2*x*chebyshev-previous
    assert direct == (x-1)**(n-1)*(n*x+(n-1))*coefficient.numerator/(n*coefficient.denominator)
    split = -Fraction(n-1,n)
    def difference(value):
        return number(coefficient)*(value-1)**(n-1)*(value-number(split))
    for index in (1,35,11431,11435):
        height = acb.zeta_zero(index).imag
        for d in (724,804):
            angle = height*(geometry["theta"](d)-geometry["theta"](244))
            delta = abs(difference((6*angle).cos()))
            assert delta>0
            if index==11431 and d==724:
                assert arb("2.6854e-12") < delta < arb("2.6855e-12")
            print("零点编号、弦类、两端面权的差",index,d,delta.str(11))
    for resolution,low,high in (("1e-6",".35825",".35827"),("1e-5",".20481",".20483")):
        cuts = []
        landmarks = [Fraction(-1),split,-Fraction(n-2,n),Fraction(1)]
        threshold = arb(resolution)
        for a,b in zip(landmarks,landmarks[1:]):
            left,right = number(a),number(b)
            def residual(z):
                return abs(difference(z))-threshold
            positive_left = residual(left)>0
            assert residual(left)*residual(right)<0
            for _ in range(48):
                middle = (left+right)/2
                value = residual(middle)
                assert value>0 or value<0
                if (value>0)==positive_left:
                    left=middle
                else:
                    right=middle
            cuts.append(left.union(right))
        a,b,c = cuts
        density = (arb.pi()-a.acos()+b.acos()-c.acos())/arb.pi()
        assert arb(low) < density < arb(high)
        print("指定面权分辨率、全部正零点中能区分两端的极限比例",resolution,density.str(12))
    # The phase-zero value is known from face normalization. The other data
    # are the actual zero probes. A real degree-18 trigonometric polynomial
    # at angular distance d from a maximum has magnitude >= max*cos(18*d)
    # for 18*d < pi/2 (the sharpened Bernstein inequality).
    speed = 6*(geometry["theta"](724)-geometry["theta"](244))
    samples = [(speed*z.imag).cos().acos() for z in acb.zeta_zeros(1,455)]
    def covering_radius(count):
        angles = [arb(0)]+sorted(samples[:count], key=float)
        assert all(a<b for a,b in zip(angles,angles[1:]))
        gaps = [(b-a)/2 for a,b in zip(angles,angles[1:])]+[arb.pi()-angles[-1]]
        return max(d.lower() for d in gaps).union(max(d.upper() for d in gaps))
    for count, amplification in ((168,"2"),(320,"1.25"),(455,"1.1")):
        threshold = (1/arb(amplification)).acos()/18
        radius = covering_radius(count)
        assert radius < threshold < covering_radius(count-1)
        print("首次通过网距判据的零点前缀、全角度误差放大上界", count, (1/(18*radius).cos()).str(10))


if __name__ == "__main__":
    verify()
