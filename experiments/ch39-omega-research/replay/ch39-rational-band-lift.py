"""有理数运算验证整条共同谱隙上的多项式次数与有界质量滤波。"""
from fractions import Fraction as F
from itertools import combinations
from math import isqrt, prod


def verify():
    # Square classes of the three algebraic points defining the two clocks.
    classes = [-d*(3468-d) for d in (244,724,804)]
    for k in range(1,4):
        for terms in combinations(classes,k):
            n = prod(terms)
            assert n < 0 or isqrt(n)**2 != n
    a,b,epsilon = F(87,10),F(48,5),F(1,10**35)
    z = 2*b/a-1
    old,new = F(1),z
    degree = 1
    while epsilon*(1+new)/2 < 84:
        old,new = new,2*z*new-old
        degree += 1
    assert degree == 137 and epsilon*(1+old)/2 < 84 <= epsilon*(1+new)/2

    # x=(20 lambda-183)/537. The two bands lie in [-1,-a0] and [a0,1].
    floor = F(3,179)
    assert (20*a-183)/537 == -floor and (20*b-183)/537 == floor
    assert (20*36-183)/537 == 1
    denominator = 10**65
    for _ in range(6):
        new_square = 4*floor/(1+floor)**2
        integer = isqrt(denominator**2*new_square.numerator//new_square.denominator)
        next_floor = F(integer,denominator)
        assert next_floor**2 <= new_square < (next_floor+F(1,denominator))**2
        floor = next_floor
    leakage = 100*(1-floor)**2/(2*(1+floor**2))
    assert leakage < F(126,10**49)
    print("纯多项式最低次数、有理步数、平方质量泄漏上界",degree,6,float(leakage))


if __name__ == "__main__":
    verify()
