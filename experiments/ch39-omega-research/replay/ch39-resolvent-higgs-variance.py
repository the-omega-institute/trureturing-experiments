"""对单个 Yukawa 顶角的四个轻双重态实分量作自由高斯平均，认证一圈真空读数。"""
from math import comb, factorial
from pathlib import Path
from runpy import run_path

from flint import arb, ctx


vacuum = run_path(str(Path(__file__).with_name("ch39-resolvent-fermion-vacuum.py")))
filtered = vacuum["filtered"]
model = filtered["model"]


def moments(h, sigma, n):
    # e^(σ²Δ/2) x_a |x|^(2n)，提出平均方向 h_a 后的有限多项式。
    return sum((comb(n,j)*factorial(n+2)//factorial(n+2-j)) * (2*sigma**2)**j * h**(2*(n-j))
               for j in range(n+1))


def average_response(h, sigma):
    moment = [moments(h,sigma,n) for n in range(6)]
    def response(left,right):
        massive = [m for m in (left,right) if m>0]
        if not massive:
            assert left.is_zero() and right.is_zero()
            return arb(1)
        n = 2*len(massive)
        divisor = arb(1)
        for m in massive:
            divisor *= m**4
        top = moment[n]/divisor
        low = top - 2*sum(1/m**2 for m in massive)*moment[n+1]/divisor
        assert low>0
        return low.union(top)
    return response


def verify():
    ctx.prec = 640
    products = filtered["products"]
    h = arb("246.22")
    scale = arb(1000)
    planck = arb("1.22089e19")
    density = arb("1e-122")*planck**4/(8*arb.pi())
    for n,w,z in ((22027,model["WITNESS"],model["HIGGS"]),
                   (28325,products["SECOND_WITNESS"],products["SECOND_HIGGS"])):
        parts,matrices = filtered["source"](n,model["column"](w),model["column"](z))
        lower, upper = (("1121.60","1121.61") if n == 22027 else ("1127.31","1127.32"))
        for width in ("0","246","1000",lower,upper,"1200"):
            sigma = arb(width)
            _,logarithm,_ = vacuum["trace_and_log"](parts,matrices,h,scale,2,average_response(h,sigma))
            ratio = logarithm/(16*arb.pi()**2*density)
            if width in ("0","246","1000",lower):
                assert ratio < 1
            else:
                assert ratio > 1
            print(n,width,"给定轻 Higgs 涨落标准差时额外费米子势负值/目标密度",ratio.str(9))


if __name__ == "__main__":
    verify()
