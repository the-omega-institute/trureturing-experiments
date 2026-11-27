"""已知一条大正权低线时，从该线的实谱权反解弱差零点的截止能标。"""
from inspect import getclosurevars
from pathlib import Path
from runpy import run_path

from flint import arb, ctx


def verify(root):
    def lib(name):
        return run_path(str(root/("ch39-"+name+".py")))

    _,spectrum = lib("adjacent-ordered-charge-six-cutoffs")["verify"](root,refine_node=2_306_888_100_086)[2]
    ctx.prec = 256
    mz,mp = arb("91.1876"),arb("1.22089e19")
    p0,p1 = (mz/mp)**2,(arb("1.176397e14")/mp)**2
    join = arb(10)**10
    fine = lib("native-rh-ray-pvm")["weak_series"](56)
    series = getclosurevars(fine).nonlocals
    coeff,ref = series["coefficients"],series["reference"]
    large = lib("full-neutral-free-cutoff")["weak_column"]
    maximum_inverse = arb(".784")
    assert p0/spectrum[0] < maximum_inverse

    def poly(inv):
        value = arb(0)
        for a in reversed(coeff):
            value = (value+a)*inv
        return value

    for i in range(8):
        cell = arb(i).union(arb(i+1))*maximum_inverse/8
        slope = arb(0)
        for k in range(len(coeff),0,-1):
            slope = slope*cell+k*coeff[k-1]
        assert slope > arb(".01567")
    a,b,c = coeff[:3]
    terms = (mz**2*a/(-ref),mz**2*b/a,mz**2*(c/a-(b/a)**2)*(-ref)/a)
    integers = (1839,17731,-283766)
    assert all(abs(v-n) < arb(".5") for n,v in zip(integers,terms,strict=True))
    def prediction(w):
        return (1839*w+17731-283766/w).sqrt()

    slack = arb(".65")
    # 弱差等于参考值加所选单线校正和其余线的加权余量，谱权之和精确为一。
    other = 1500*(large(arb("1e-12"),p0,p1)-ref).abs_upper()
    series_error = series["remainder"] + 2*join*series["p_ratio"]
    assert ref+500*(large(join*p0,p0,p1)-ref).abs_upper()+other < 0
    domain = [(arb(100),arb(500))]
    cells = 0
    while domain:
        left,right = domain.pop()
        # 预测能量随 w 递增；P(s) 随 s 递增。两个端点分别比较整个权重区间。
        qlo,qhi = prediction(right)-slack,prediction(left)+slack
        s_lo,s_hi = (mz/qlo)**2,(mz/qhi)**2
        assert 1/join < s_lo < maximum_inverse and 1/join < s_hi < maximum_inverse
        low_corr,high_corr = poly(s_lo)-series_error,poly(s_hi)+series_error
        assert low_corr > 0 and high_corr > 0
        if ref+left*low_corr-other > 0 and ref+right*high_corr+other < 0:
            cells += 1
            continue
        assert right-left > arb("1e-5")
        mid = ((left+right)/2).mid()
        domain.extend(((left,mid),(mid,right)))
    print("系数的原值、向最近整数取整及覆盖闭段数",[v.str(14) for v in terms],integers,cells,flush=True)
    for j,i in enumerate((0,15,0)):
        weights = next(r[2][2] for r in lib("native-rh-ray-pvm")["shell_weights"](root,2_306_888_100_092+j) if r[0]==6)
        w = weights[i].real
        assert 100 < w < 500 and sum(abs(v.real) for k,v in enumerate(weights) if k!=i)<1500
        print("壳、谱线、实谱权与有理式预测的 GeV 能量",j+7,i+1,w.str(13),prediction(w).str(14),flush=True)


if __name__ == "__main__":
    verify(Path(__file__).resolve().parent)
