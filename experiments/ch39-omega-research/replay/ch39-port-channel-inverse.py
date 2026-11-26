"""Recover the channel count, excited energy and contact from interior clock moments."""
from fractions import Fraction as Q

from flint import arb, ctx


def verify():
    ctx.prec = 144
    count = 624
    weight = Q(count,count+1)
    mean = 2+weight
    second = 6+5*weight
    third = 20+25*weight
    raw1 = mean-2
    raw2 = second-6-4*raw1
    raw3 = third-20-18*raw1-6*raw2
    assert raw1*raw3 == raw2**2
    assert raw1**2/(raw2-raw1**2) == count and raw2/raw1 == 1
    at_contact = 1+weight+1
    assert at_contact-mean+1 == 1

    def real(q):
        return arb(q.numerator)/q.denominator
    def distance(other):
        # Nearest point on the competing parabola lies on a square corner.
        direction = 1 if other > count else -1
        def signed_gap(radius):
            x = real(mean)+direction*radius
            y = x**2+2+(x-2)**2/other
            return direction*(y-real(second))+radius
        lo, hi = arb(0),arb(".000001")
        assert signed_gap(lo) < 0 < signed_gap(hi)
        for _ in range(44):
            mid = ((lo+hi)/2).mid()
            gap = signed_gap(mid)
            if gap < 0:
                lo = mid
            else:
                assert gap > 0
                hi = mid
        return lo.union(hi)
    lower, upper = distance(count-1), distance(count+1)
    threshold = upper/2
    assert upper < lower
    assert arb("1.825e-7") < threshold < arb("1.83e-7")
    print("native first three moments", mean,second,third,"and contact mean",at_contact)
    print("nearest other count 625, sharp common per-reading error threshold",threshold.str(14))


if __name__ == "__main__":
    verify()
