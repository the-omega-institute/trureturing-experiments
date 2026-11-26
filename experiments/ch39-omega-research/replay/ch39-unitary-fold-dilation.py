"""从严格九十度折叠压缩出三代对角 Yukawa 的最少空间维数。"""
from fractions import Fraction
from math import isqrt, prod
from pathlib import Path
from runpy import run_path

from sympy import factorint

gram = run_path(str(Path(__file__).with_name("ch39-native-shell-rational-gram.py")))


def minimum_dimensions(left_entries, right_dim):
    # Values are the normalized diagonal singular values, with extra zero
    # columns when the original left space is larger than the right space.
    defect = [1-v*v for v in left_entries if v!=1]
    rank = len(defect)
    determinant = prod(defect, start=Fraction(1))
    primes = {2}
    for value in defect:
        primes.update(factorint(value.numerator))
        primes.update(factorint(value.denominator))
    bad, two_blocked, inert = [],[],[]
    for p in sorted(primes):
        exponent, unit = gram["local_unit"](-determinant,p)
        if p%4==3 and exponent%2:
            inert.append(p)
        parity = prod(gram["hilbert"](x,y,p) for j,x in enumerate(defect) for y in defect[:j])
        target = parity*gram["hilbert"](determinant,-1,p)
        if target == -1:
            bad.append(p)
            square_unit = (unit==1) if p==2 else (pow(unit,(p-1)//2,p)==1)
            if exponent%2==0 and square_unit:
                two_blocked.append(p)
    rational_square = all(isqrt(v)**2==v for v in (determinant.numerator,determinant.denominator))
    if two_blocked:
        rational_extra = 3
    elif bad:
        rational_extra = 2
    else:
        rational_extra = int(not rational_square)
    base = right_dim+rank
    return (base,base+int(bool(inert)),base+rational_extra),two_blocked


def verify():
    # Treat just the printed central values as diagonal rational test data.
    data = list(map(Fraction, ".00000615 .00311 .8616 .0000135 .000268 .01401 .0000028227 .00059465 .010101".split()))
    up,down,charged = [[v/data[2] for v in data[i:i+3]] for i in (0,3,6)]
    neutrino = [Fraction(0)]*3
    all_values = 3*up+3*down+charged+neutrino
    minimal,obstruction = minimum_dimensions(all_values,21)
    assert minimal == (42,43,45) and obstruction == [7]
    blocks = [minimum_dimensions(values,right)[0]
              for values,right in ((up,3),(down,3),(charged,3),(neutrino,0))]
    assert blocks == [(5,6,7),(6,7,8),(6,7,8),(3,3,3)]
    split = [3*(a+b)+c+d for a,b,c,d in zip(*blocks)]
    complete = [6*max(a,b)+2*max(c,d) for a,b,c,d in zip(*blocks)]
    assert split==[42,49,56] and complete==[48,56,64]
    print("每侧最少方向（实数、Gaussian 有理数、有理数）：总压缩、分电荷色、左弱双重态配齐",
          minimal,split,complete)
    # An exact rational rotation for each mass when the printed experimental
    # intervals are used. Choose the circle parameter by integer square root.
    denominator = 10**12
    errors = []
    for value in up+down+charged:
        radicand = (1-value)/(1+value)
        numerator = isqrt(denominator**2*radicand.numerator//radicand.denominator)
        parameter = Fraction(numerator,denominator)
        cosine = (1-parameter**2)/(1+parameter**2)
        sine = 2*parameter/(1+parameter**2)
        assert cosine*cosine+sine*sine == 1
        errors.append(abs(value-cosine)*data[2])
    assert max(errors) < Fraction(2,10**12)
    print("有限观测精度下用有理二阶旋转实现的最大绝对 Yukawa 误差", float(max(errors)))


if __name__ == "__main__":
    verify()
