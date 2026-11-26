"""高处标准逆平方尾权在任意固定相位接受窗内的两个归一时间读数。"""
from pathlib import Path
from runpy import run_path

from flint import arb, ctx


clock = run_path(str(Path(__file__).with_name("ch39-arithmetic-height-bandwidth.py")))


def uniform_reading(argument):
    # 对有限有理被积式在 [0,1] 逐项积分，奇偶截断给出两侧。
    square = argument**2
    even = sum((-square)**j/(2*j+1) for j in range(17))
    odd = even-square**17/35
    return odd.union(even)


def verify():
    ctx.prec,ctx.cap = 1536,3
    target = clock["times"]()
    ctx.prec = 192
    cases = (("首项校准",".06943943437405956",".06943943437405958"),
             ("两项最大误差最小",".06935860590438615",".06935860590438618"))
    for name,low,high in cases:
        lo,hi = arb(low),arb(high)
        for endpoint,sign in ((lo,1),(hi,-1)):
            residual = uniform_reading(endpoint)-target[1]
            if name=="两项最大误差最小":
                residual += uniform_reading(2*endpoint)-target[2]
            assert sign*residual > 0
        scale = lo.union(hi)
        first = uniform_reading(scale)-target[1]
        second = uniform_reading(2*scale)-target[2]
        if name=="首项校准":
            assert arb("1.83377e-5") < -second < arb("1.83378e-5")
        else:
            assert arb("3.71810e-6") < first < arb("3.71812e-6")
            assert arb("3.71810e-6") < -second < arb("3.71812e-6")
        print(name,"尺度与第二项误差",scale.str(17),second.str(12))


if __name__=="__main__":
    verify()
