"""在有界隐标量势下，检验带电费米子平方 Yukawa 与中性标量门户的跑动。"""
from fractions import Fraction
from itertools import permutations, product
from pathlib import Path
from runpy import run_path

from flint import arb, ctx


sm = run_path(str(Path(__file__).with_name("ch39-two-loop-shell-separation.py")))


def finite_moment_check():
    # 两坐标例子检验有限原子恒等式。坐标值 0,±3 给方差 3，再按 1,2 缩放。
    atoms = ((0,Fraction(2,3)),(-3,Fraction(1,6)),(3,Fraction(1,6)))
    for i in range(5):
        total = sum(pa*pb*a**i*(2*b)**(4-i) for (a,pa),(b,pb) in product(atoms, repeat=2))
        assert total == (432,0,36,0,27)[i]


def verify():
    ctx.prec = 192
    finite_moment_check()
    rel = Fraction(1,1000)
    g1, _, g3 = sm["hull_inputs"](".467453 .63811 1.0583",rel)
    ys = sm["hull_inputs"](".00000615 .00311 .8616 .0000135 .000268 .01401 .0000028227 .00059465 .010101",rel)
    p = [y**2 for y in ys]
    hyper = 3*g1**2/5
    gauge = [16*g3**2+4*hyper/3]*3 + [16*g3**2-2*hyper/3]*3 + [6*hyper]*3
    multiplicities = [12]*6+[4]*3
    weights = [d*v**2 for d,v in zip(multiplicities,p)]
    total = sum(weights)
    independent = sum(w*(v+c) for w,v,c in zip(weights,p,gauge))
    costs = []
    for permutation in permutations(range(3)):
        costs.append(3*sum(weights[i]*p[3+j]+weights[3+j]*p[i] for i,j in enumerate(permutation)))
    minimum = min(v.lower() for v in costs)
    bound = 2*(independent+minimum)/(16*arb.pi()**2*total)
    assert bound > arb(".236")
    # 隐四次与 λ 均可取零，反序置换给定 Yukawa 时确实最小。
    assert p[0] < p[1] < p[2] and p[3] < p[4] < p[5]
    # 对同一组输入，交换逆序两列的代价差是两项正数之和。不同区间无需彼此分离。
    print("平方门户迹与费米子四次矩之比的对数斜率下界",bound.lower().str(11))


if __name__=="__main__":
    verify()
