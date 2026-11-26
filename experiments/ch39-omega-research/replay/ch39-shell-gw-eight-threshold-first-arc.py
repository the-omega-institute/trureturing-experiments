"""First continuous returns to either of the two eighthreshold cones of the shell clocks."""
from math import gcd

from flint import acb, arb, ctx


def main():
    ctx.prec = 320
    q = arb
    angles = [(1-q(j)/26).acos() for j in (3,5,7)]
    f,g = [a-angles[0] for a in angles[1:]]
    ratio, w = g/f, q("1e-18")
    k,m = 640942321255267281,1173203932801135019
    assert gcd(k,m) == 1 and k%2 == m%2 == 1
    residual = k*ratio-m
    assert abs(residual) < q("3.329e-19")
    # Use integer multiples of pi for the union of the even-even and odd-odd centers.
    tol = (1+ratio)*w/q.pi()
    assert k*(abs(residual)+tol) < q(".792")
    assert 0 < w/g < 1 < acb.zeta_zero(1).imag
    for parity in (1,2):
        step = parity*q.pi()
        offset = step*residual
        # u = first angular deviation; second deviation = ratio*u + offset.
        enter = (-w-offset)/ratio
        leave = -offset/(1+ratio)
        assert 0 < enter < leave < w
        assert ratio*w+offset < w
        # Positive v >= |u| is impossible for all u <= w once (ratio-1)w+offset<0.
        assert (ratio-1)*w+offset < 0
        for t in (enter,leave):
            assert abs(ratio*t+offset) > q("1e-19")
        h = (step*k+enter)/f
        duration = (leave-enter)/f
        if parity == 1:
            assert q("1.3844650454655219e19") < h < q("1.3844650454655221e19")
            assert q("2.36e-18") < duration < q("2.38e-18")
        else:
            assert q("2.7689300909310438e19") < h < q("2.7689300909310440e19")
            assert q("9.81e-19") < duration < q("9.82e-19")
        print("first return to", "pi,pi" if parity == 1 else "0,0",
              "height",h.str(55),"interval length",duration.str(15))


if __name__ == "__main__":
    main()
