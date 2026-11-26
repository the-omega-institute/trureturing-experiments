"""三轻模与十六重模的壳热迹能在九枚实际零点上产生多大的相邻谱权增长。"""
from pathlib import Path
from runpy import run_path

from flint import acb, arb, ctx


root = Path(__file__).parent
weights = run_path(str(root / "ch39-arithmetic-weight-rise.py"))
shell = run_path(str(root / "ch39-three-ray-chiral-boundary.py"))


def heavy_spectrum(number):
    h = acb.zeta_zero(number).imag
    neutral = shell["laplacian"](0,h)
    center = sum(neutral[i,i].real for i in range(19))/18
    bands = {}
    for q in (1,3,2,4,6):
        eigen = shell["laplacian"](q,h).eig(multiple=True)
        levels = sorted((abs(z.real*(z.real-center)) for z in eigen),key=float)
        assert levels[2] < 9 < levels[3]
        bands[q] = levels[3:]
    return bands


def partition(bands, t):
    return 3+sum((-t*b).exp() for b in bands)


def below_ratio(old, new, ratio):
    # 全区间对 log m 的 Lipschitz 界常常已足够。
    crude = sum(max(arb(0), (a/b).log().upper()) for a,b in zip(old,new)) / arb(1).exp()
    if 3*(ratio-1) > crude:
        return
    unit = min(new).lower()
    a,b = [[v/unit for v in row] for row in (old,new)]
    high = arb(8)
    assert 16*(-high).exp() < 3*(ratio-1)
    def verify(lo, hi, depth):
        mid, radius = (lo+hi)/2, (hi-lo)/2
        exp_a = [(-mid*x).exp() for x in a]
        exp_b = [(-mid*x).exp() for x in b]
        difference = 3+sum(exp_b) - ratio*(3+sum(exp_a))
        slope = ratio*sum(x*v for x,v in zip(a,exp_a)) - sum(x*v for x,v in zip(b,exp_b))
        second = sum(x*x*(-lo*x).exp() for x in b)
        second += ratio*sum(x*x*(-lo*x).exp() for x in a)
        if difference + abs(slope)*radius + second*radius**2/2 < 0:
            return
        assert depth < 30
        verify(lo,mid,depth+1)
        verify(mid,hi,depth+1)
    verify(arb(0),high,0)


def verify():
    ctx.prec = 224
    spectra = [heavy_spectrum(number) for number in weights["numbers"]]
    for power,bound,temp,lower in ((1,"1.007",".0522","1.00694"), (2,"1.01148",".0012","1.01145")):
        for old,new in zip(spectra,spectra[1:]):
            for q in old:
                below_ratio([v**power for v in old[q]], [v**power for v in new[q]], arb(bound))
        a,b = [[v**power for v in row[3]] for row in spectra[-2:]]
        seen = partition(b,arb(temp))/partition(a,arb(temp))
        assert seen > arb(lower)
        print("质量幂", power, "相邻最大倍率下界样例与全域上界",seen.str(10),bound)
    assert arb("1.01148")**12 < arb("1.14827")


if __name__ == "__main__":
    verify()
