"""固定谱线分布，允许目标独自改变标尺时，两引力矩的迁移障碍。"""
from itertools import permutations
from math import factorial
from pathlib import Path
from runpy import run_path

from flint import arb, arb_mat

ROOT = Path(__file__).resolve().parent


def derivative(coeff, s):
    """在原标尺的 r 基底中表示对目标标尺比 s 的导数。"""
    numbers = coeff.entries()
    terms = [x*(2-i)*s**(1-i) for i,x in enumerate(numbers[:-1])]
    terms[2] = numbers[-1]/s
    return arb_mat([[*terms, 0]])


def verify(*, forward=None):
    load = run_path(str(ROOT/"ch39-four-background-crossing.py"))["load"]
    order_check = load("fifth-shell-vacuum-order")
    old = load("fifth-shell-frozen-witnesses")["NUMBERS"]
    numbers = (*old, *(order_check["HOLDOUT"] if forward is None else forward))
    ctx = load("shared-large-scale-gap")["load"]()
    fields, kernel = ctx[1], ctx[-2]
    scale = arb(7)/1000
    shell = load("fourth-fifth-shell-exchange")["FIFTH_SHELL"]
    edges = load("shared-twelve-node-gauge-gap")["paired_lines"](
        ROOT, load("native-rh-ray-pvm")["shell_weights"], shell)
    moments = load("fifth-shell-transfer")["shell_moments"](edges, kernel, scale)
    _, lib, spectra = load("two-zero-cutoff-gravity")["gravity_fields"](ROOT,numbers,tight_neutral=True)
    assert all(max(z.upper() for _,_,z in rows) < min(arb(1),scale*min(a.lower() for a,_,_ in edges)/5)
               for rows in spectra)
    assert all([g for g,_,_ in rows]==[g for g,_,_ in spectra[0]] for rows in spectra)
    series = [kernel["series"](rows,fields) for rows in spectra]
    newton = [arb_mat([[1,0,0,0,0]])*mat for mat in series]
    vacuum = [arb_mat([[0,1,0,0,0]])*mat for mat in series]
    weight = sum(abs(w) for _,_,w in edges)
    n = kernel["ORDER"]+1
    # z/(sx)<2, s>=1/10。尾的第 j 项对 s 求导后乘 j/s；
    # 核的第 j 项均含 1/j，导数余项因子至多 10。
    tail = 2*arb(2)**n*weight/(4*arb.pi()*factorial(n))
    kt = [tail*sum(abs(2*a) for _,a,_ in rows) for rows in spectra]
    wt = [tail*sum(abs(g) for g,_,_ in rows) for rows in spectra]
    eps = arb("1e-122")
    relative_small = arb(1)/10
    def massless_error(index,s):
        rows = spectra[index]
        # x(1-E2(z/x)) 随 x 增，导数为 1-exp(-z/x)。
        def at(right):
            x = scale*s*right
            return x*sum(abs(a)*(1-lib["expint"](2,z/x)) for _,a,z in rows)
        return sum(abs(w)*at(b) for _,b,w in edges)/(2*arb.pi())
    tiny = [massless_error(i,relative_small) for i in range(len(numbers))]
    source = [massless_error(i,1) for i in range(len(old))]
    def expect(row):
        return (row*moments)[0,0]
    pairs = [(3,j) for j in range(4,len(numbers))]
    if forward is None:
        pairs = [*permutations(range(len(old)),2), *pairs]
    corrections, gaps = [], []
    for i,j in pairs:
        small = relative_small*(1+source[i])+tiny[j]
        assert small < 1
        slopes = []
        left = 1
        while scale*left/10 < 1:
            s = (arb(left)/10).union(arb(2*left)/10)
            slopes.append((1+expect(derivative(newton[j],s)-newton[i])-kt[i]-10*kt[j]).lower())
            left *= 2
        slope = min(slopes)
        assert slope > arb(".98")
        kdelta = abs(expect(newton[j]-newton[i]))+kt[i]+kt[j]
        assert 1-kdelta+slope > 1 and 2*scale < 1
        step = (kdelta/slope).upper()
        assert step < arb(".001")
        s = arb(1)+arb(0, step)
        wd = derivative(vacuum[j],s)-2*s*vacuum[i]
        # 二次项系数仅依赖总 g，两处完全相等。
        wd[0,0] = 0
        multiplier = wd[0,1]/newton[i][0,1]
        err = 10*wt[j]+2*abs(s)*wt[i]+abs(multiplier)*kt[i]
        speed = abs(2*s*eps+multiplier+expect(wd-multiplier*newton[i]))+err
        # 独立重算固定标尺的真空差，消去 x 项。
        dw = vacuum[j]-vacuum[i]
        ratio = dw[0,1]/newton[i][0,1]
        at_one = ratio+expect(dw-ratio*newton[i])
        margin = abs(at_one)-wt[i]-wt[j]-abs(ratio)*kt[i]-speed*step
        print(numbers[i],numbers[j],"小标尺牛顿上界、最大相对调整、调整后真空差下界",
              small.upper().str(7),step.str(8),margin.lower().str(9),flush=True)
        assert margin > 0
        corrections.append(step)
        gaps.append(margin.lower())
    if forward is None:
        assert max(corrections)<arb(".000748") and min(gaps)>arb("2.416e-9")
    if forward==tuple(range(2345724,2345747,2)):
        assert min(gaps)>arb("3.0886e-8")
    print("全部配对最大相对调整、最小真空差下界",max(corrections).str(9),min(gaps).str(9))


if __name__=="__main__":
    verify()
