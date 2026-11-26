"""Exact finite-precision clock lifts on the four-dimensional shell representation."""
from fractions import Fraction as Q
import json
from pathlib import Path


def turn(chart, den, scale):
    sign, k = chart
    t = Q(k,den)
    q = 17*scale**2
    return sign*(1-q*t*t)/(1+q*t*t), 2*sign*scale*t/(1+q*t*t)


def reference(chart, den):
    sign, k = chart
    t = Q(k,den)
    return sign*(1-t*t)/(1+t*t), 2*sign*t/(1+t*t)


def readings(row, den, scale):
    _, _, _, a, b = row
    x, u = turn(a,den,scale)
    y, v = turn(b,den,scale)
    assert x*x + 17*u*u == y*y + 17*v*v == 1
    return x,y,x*y+17*u*v,x*y-17*u*v


def targets(row, den):
    _, _, _, a, b = row
    x,u = reference(a,den)
    y,v = reference(b,den)
    return x,y,x*y+u*v,x*y-u*v


def main():
    data = json.loads(Path(__file__).with_name("ch39-five-tick-record.json").read_text())
    den = 2**data["bits"]
    rows = data["rows"]
    left = Q(43581196642368,179689877047297)
    right = Q(223271073689665,920570219967553)
    scale = Q(left.numerator+right.numerator,left.denominator+right.denominator)
    assert right.numerator*left.denominator-left.numerator*right.denominator == 1
    worst = max(abs(x-y) for r in rows for x,y in zip(readings(r,den,scale),targets(r,den)))
    tolerance = Q(1,2**96)
    assert Q("0.6520")*tolerance < worst < Q("0.6521")*tolerance

    # One single clock forces q=17*scale^2 into (1/2,2), since the
    # sign-corrected real part is strictly decreasing for nonzero chart t.
    row = rows[40]
    t,u = Q(row[3][1],den),Q(row[4][1],den)
    real = lambda q: (1-q*t*t)/(1+q*t*t)
    assert real(Q(1,2)) > real(1)+tolerance and real(2) < real(1)-tolerance
    # On that interval q/[(1+q*t^2)(1+q*u^2)] is increasing:
    # its divided difference has the sign of 1-q1*q2*t^2*u^2.
    # The relative phase real part is sign_a*sign_b times
    # 1 - 2*(t-u)^2*q/[(1+q*t^2)(1+q*u^2)].
    assert 2*abs(t*u) < 1 and row[3][0]*row[4][0] == -1 and t != u
    exact = targets(row,den)[2]
    assert Q(1,2) < 17*left**2 < 17*right**2 < 2
    assert readings(row,den,left)[2] < exact-tolerance
    assert readings(row,den,right)[2] > exact+tolerance
    assert scale.denominator == 1100260097014850

    # Opposite ways to lift both recorded clocks fit every even reading.
    # The two signed loop averages differ by -(14/68)*sum(weight*b_j).
    total = sum((r[2]*turn(r[3],den,scale)[1]) // 1 for r in rows)
    assert Q(".02065424526") < -Q(14*(total+46),68*den)
    assert -Q(14*total,68*den) < Q(".02065424527")
    print("Four rational dimensions. Minimizing denominator for this common positive Cayley scale:", scale)
    print("184 real phase errors below 0.6521 times 2^-96. Opposite signed loop readings differ by 0.02065424526...")


if __name__ == "__main__":
    main()

