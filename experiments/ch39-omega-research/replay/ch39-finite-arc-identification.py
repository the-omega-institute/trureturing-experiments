"""两枚相邻编号零点在有界整壳上的弧差分离与无界的非十七倍逼近。"""
from fractions import Fraction
from math import isqrt, lcm

from flint import acb, arb, ctx, fmpz_poly


N = 2_306_888_100_085
BASE = 867
REF = 745
ENDS = (505, 465)


def ar(f):
    return arb(f.numerator)/f.denominator


def polynomial(a, b, m):
    return fmpz_poly([a*a+b*b-m*m, -2*a*b, m*m])


def lift(seed, k):
    # 球半径 sqrt(3)，从有理基点 (1,1,1) 取第二个直线交点。
    v = [17*k*(s-17)+1 for s in seed]
    den = sum(z*z for z in v)
    assert den % 17 == 3
    w = tuple(1-Fraction(2*sum(v)*z,den) for z in v)
    assert sum(z*z for z in w) == 3
    return w


def verify():
    ctx.prec = 224
    h0,h1 = [acb.zeta_zero(j).imag for j in (N,N+1)]
    gap = h1-h0
    assert arb(".07446577654297158521") < gap < arb(".07446577654297158522")
    assert gap*arb.pi() < arb.pi()
    denom = 8*81*BASE**4
    assert denom == 366143556563208
    assert gap/denom > arb("2.0337e-16")
    assert 2*gap/(arb.pi()*denom) > arb("1.2947e-16")
    assert arb("13.42898773670795") < 1/gap < arb("13.42898773670796")
    print("编号、高度差、每弧度相位增量的角频误差倍率", N, gap.str(22), (1/gap).str(19))
    target = []
    for b in ENDS:
        P = polynomial(REF, b, BASE)
        discriminant = int(P[1]**2-4*P[2]*P[0])
        assert int(P.content()) == 1 and isqrt(discriminant)**2 != discriminant
        alpha = (arb(b)/BASE).acos()-(arb(REF)/BASE).acos()
        assert arb(".4") < alpha < arb(".5") and P(alpha.cos()).contains(0)
        A = acb(0,alpha*h0).exp()
        B = acb(0,alpha*h1).exp()
        assert (B/A-acb(0,alpha*gap).exp()).contains(0)
        target.append(alpha)
        print("原钟端、最小多项式系数、相位增量角",b,list(P),(gap*alpha).str(20))
    for k in (1, 100, 10000):
        outer = lift((1,5,29),k)
        near = lift((13,13,23),k)
        # 交换最后两坐标将原近点映到第二端，交换第一和第三坐标映到第三端。
        near_points = (near, (near[0],near[2],near[1]), (near[2],near[1],near[0]))
        points = (outer,*near_points)
        u = lcm(*(v.denominator for p in points for v in p))
        assert u % 17 != 0
        integer_points = [tuple(int(u*v) for v in p) for p in points]
        assert all(min(p) >= 0 and sum(t*t for t in p) == 3*u*u for p in integer_points)
        cosines = [sum(a*b for a,b in zip(outer,v))/3 for v in near_points]
        frequencies = [ar(c).acos()-ar(cosines[0]).acos() for c in cosines[1:]]
        for ref, c, frequency, original, end in zip([cosines[0]]*2,cosines[1:],frequencies,target,ENDS):
            m = 3*u*u
            coefficients = (int(ref*m),int(c*m))
            assert all(v.denominator == 1 for v in (ref*m,c*m))
            P, Q = polynomial(REF,end,BASE),polynomial(*coefficients,m)
            assert P.resultant(Q) != 0 and not (frequency-original).contains(0)
        if k == 10000:
            errors = [abs(gap*(v-a)) for v,a in zip(frequencies,target)]
            assert errors[0] < arb("2.785e-7") and errors[1] < arb("3.102e-7")
            assert u == 10298688979560253193880003
        print("近似参数、共同整壳 u、两枚增量角的两钟误差",
              k, u, *[(gap*(v-a)).str(13) for v,a in zip(frequencies,target)])


if __name__ == "__main__":
    verify()
