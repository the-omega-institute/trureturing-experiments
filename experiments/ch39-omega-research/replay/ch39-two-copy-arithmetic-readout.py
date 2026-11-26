"""两份坐标轮换不变的轻态无误识别九枚实际零点的最优共同效率。"""
from pathlib import Path
from runpy import run_path

from flint import acb_mat, arb, ctx


root = Path(__file__).parent
mix = run_path(str(root / "ch39-two-boundary-site-mixing.py"))
numbers = (8359,11661,11934,12207,22027,28325,48067,58214,95478)


def exceeds(overlaps, copies, probability):
    # G_k - pI 的首个严格负主子式排除半正定；顺序主子式全正保证正定。
    matrix = acb_mat([[1-probability if i == j else overlaps[i,j]**copies
                       for j in range(9)] for i in range(9)])
    for n in range(2,10):
        det = acb_mat([[matrix[i,j] for j in range(n)] for i in range(n)]).det().real
        if det < 0:
            return False
        assert det > 0
    return True


def verify():
    ctx.prec = 512
    frames = [mix["light_modes"](n) for n in numbers]
    # 荷一的轮换不变轻线，整列的任意相位不影响 Gram 谱。
    states = acb_mat([[f[i,0] for f in frames] for i in range(19)])
    adjoint = mix["base"]["adjoint"]
    overlaps = adjoint(states)*states
    gram = acb_mat([[overlaps[i,j]**2 for j in range(9)] for i in range(9)])
    inverse = gram.inv()
    identity = acb_mat([[int(i==j) for j in range(9)] for i in range(9)])
    assert all(z.contains(0) for z in (gram*inverse-identity).entries())
    scaled = inverse * arb("1e-21")
    power = scaled**8
    trace8 = sum(power[i,i].real for i in range(9))
    power = power*scaled
    trace9 = sum(power[i,i].real for i in range(9))
    lo = arb("1e-21")/trace8.root(8)
    hi = arb("1e-21")*trace8/trace9
    assert arb("2.37133020972e-22") < lo < hi < arb("2.37133020974e-22")
    print("两份轮换不变轻态的最优共同无误识别概率",lo.str(14),hi.str(14))
    for probability, first in ((arb(1)/2,19458133657),(arb(9)/10,64372440933)):
        assert not exceeds(overlaps,first-1,probability)
        assert exceeds(overlaps,first,probability)
        print("共同效率达到",probability.str(3),"的最小副本数",first)



if __name__ == "__main__":
    verify()
