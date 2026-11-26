"""放开整数占据后逐面取最低实作用，在有限形变上重放 Regge 项与角亏平方修正。"""
from flint import arb, ctx


def cell(n, t, b):
    result = arb(0)
    for count, x, deficit in [(30, b**2*t*(4-t)/16, 2*arb.pi()-4*((3*t-2)/(2+2*t)).acos()),
                              (20, b**2*(4*t-1)/16, 2*arb.pi()-6*((3-2*t)/(2+2*t)).acos())]:
        h = deficit/(8*arb.pi())
        # F''>0 above 3√x/4. At all smaller p the energy exceeds that of a rounded √x.
        radius = x.sqrt()
        lo = int((n*radius*3/4).floor().unique_fmpz()) + 1
        sample = (n*radius+arb(1)/2).floor()/n
        assert (7*x/16)**2/2 - abs(h) > sample*h+(sample**2-x)**2/2
        hi = n
        while lo < hi:
            k = (lo+hi)//2
            p = arb(k)/n
            step = arb(1)/n
            # n times the adjacent-sector energy difference, avoiding cancellation.
            difference = h + (2*p+step)*(p**2+(p+step)**2-2*x)/2
            if difference > 0:
                hi = k
            else:
                assert difference < 0
                lo = k+1
        p = arb(lo)/n
        result += count*(h*p+(p**2-x)**2/2)
    return result


def verify():
    ctx.prec = 192
    q = 2*arb(3).sqrt()-3
    b = 4*q/(arb(5)/3).sqrt()
    flat = arb(2)/3
    eps = arb(1)/1000
    regge = -54*b/(arb(5).sqrt()*8*arb.pi())
    curvature_term = (30*(arb(18)/5)**2/(1-q)**2 + 20*(18*arb(3).sqrt()/5)**2/q**2)
    expected = regge-curvature_term/(4*(8*arb.pi())**2)
    print("实作用的连续占据形状 Hessian 与其中的 Regge 项",expected.str(12),regge.str(12))
    for n in (5432, 75658, 1053800):
        second = (cell(n, flat+eps,b)+cell(n,flat-eps,b)-2*cell(n,flat,b))/eps**2
        print("槽位、步长 1/1000 的实际整数中心差分",n,second.str(12))
        assert abs(second-expected) < 75/(n*eps)**2 + arb(".001")


if __name__ == "__main__":
    verify()
