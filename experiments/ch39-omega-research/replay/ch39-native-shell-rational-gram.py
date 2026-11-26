"""单位十九点壳的有理 Gram 最小嵌入维数，按局部希尔伯特符号判定。"""
from fractions import Fraction
from math import isqrt, prod
from pathlib import Path
from runpy import run_path

from flint import fmpz_mat
from sympy import factorint


def split(value,p):
    exponent = 0
    while value % p == 0:
        value //= p
        exponent += 1
    return exponent,value


def local_unit(value,p):
    n,u = split(value.numerator,p)
    m,v = split(value.denominator,p)
    modulus = 8 if p==2 else p
    return n-m,(u*pow(v,-1,modulus)) % modulus


def hilbert(a,b,p):
    n,u = local_unit(a,p)
    m,v = local_unit(b,p)
    if p==2:
        exponent = ((u-1)//2)*((v-1)//2)+n*((v*v-1)//8)+m*((u*u-1)//8)
    else:
        exponent = (n*m*((p-1)//2)+m*int(pow(u,(p-1)//2,p)!=1)
                    +n*int(pow(v,(p-1)//2,p)!=1))
    return (-1)**(exponent % 2)


def verify():
    geometry = run_path(str(Path(__file__).with_name("ch39-arc-relative-determinant.py")))
    rows = [[0]*19 for _ in range(19)]
    for i,j,_,_ in geometry["edges"]:
        rows[i][i] += 1
        rows[j][j] += 1
        rows[i][j] -= 1
        rows[j][i] -= 1
    dets = [1]+[int(fmpz_mat([row[:k] for row in rows[:k]]).det()) for k in range(1,19)]
    determinant = dets[-1]
    assert determinant == 41501049051635045807 and isqrt(determinant)**2 != determinant
    primes = sorted({2} | {int(p) for value in dets for p in factorint(value)})
    # The leading minors give a rational diagonalization by successive pivots.
    diagonal = [Fraction(b,a) for a,b in zip(dets,dets[1:])]
    def target(p):
        h = prod(hilbert(x,y,p) for i,x in enumerate(diagonal) for y in diagonal[:i])
        return h*hilbert(determinant,-1,p)
    bad = [p for p in primes if target(p)==-1]
    def square_at_p(a,p):
        n,u = local_unit(a,p)
        return n%2==0 and (u==1 if p==2 else pow(u,(p-1)//2,p)==1)
    obstruct = [p for p in bad if square_at_p(-determinant,p)]
    assert bad == [2,11,17,71] and obstruct == [2]
    print("Bad places for one and two complementary coordinates",bad,obstruct)
    print("Rational Gram dimension",21 if obstruct else 20 if bad else 19)


if __name__ == "__main__":
    verify()
