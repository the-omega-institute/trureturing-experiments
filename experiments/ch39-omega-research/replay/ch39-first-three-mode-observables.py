"""按编号取前 8359 枚实际零点上的全部五荷三轻模候选。"""
from concurrent.futures import ProcessPoolExecutor
from contextlib import redirect_stdout
from functools import cache
from io import StringIO
from pathlib import Path
from runpy import run_path
import traceback

from flint import acb, arb, ctx

root = Path(__file__).parent
joint = run_path(str(root/"ch39-chiral-joint-vacuum-density.py"))
mask = run_path(str(root/"ch39-three-mode-arithmetic-mass.py"))


@cache
def model_inputs():
    model,_,_,initial,_ = joint["chiral"]["source"]["mean"]["inputs"]()
    return model,initial


def observables(n):
    model,initial = model_inputs()
    ctx.prec=224
    levels = model["shell"]["shell_spectra"]((0,*model["charges"]),zero_index=n)
    center = sum(levels[0])/18
    layers = {q:joint["chiral"]["ordered"]([v*abs(v-center) for v in eig])
              for q,eig in levels.items()}
    return layer_observables(model,initial,layers)


def layer_observables(model,initial,layers):
    assert all(layers[q][2]<9<layers[q][3] for q in model["charges"])
    vacuum = joint["vacuum"]
    mass = 32*arb("2.435e18")/724**2*layers[2][15]
    a = (arb("246.22")**2/(mass*arb("1.22089e19")))**2
    balanced = vacuum["paired_spectrum"](vacuum["minimal"]["fields"](layers),min_gap="1e-18")
    rows = balanced+vacuum["low_sector"](4,a)
    assert vacuum["newton"](rows,arb(".03"))<1<vacuum["newton"](rows,arb(".033"))
    coefficient=(9*arb(3).log()-12*arb(2).log())/2
    curve=a**2*coefficient/(arb.pi()*arb("1e-122"))*(1-arb(0).union(a/(arb(".03")*coefficient)))
    with redirect_stdout(StringIO()):
        energy,_,alpha=joint["gauge"]["crossing_enclosure"](
            model,initial,layers,False,724,"29.04",".1184",recenter=True)
    return energy,alpha,curve


def worker(n):
    try:
        return [v.str(28,more=True) for v in observables(n)]
    except Exception:
        raise RuntimeError(f"零点 {n}\n{traceback.format_exc()}") from None


def scan(numbers):
    with ProcessPoolExecutor(max_workers=4) as pool:
        return [(n,*(arb(v) for v in values)) for n,values in zip(numbers,pool.map(worker,numbers),strict=True)]


def candidate_numbers():
    heights=mask["verify"](8359)
    numbers=[]
    for h in heights:
        left=(h-arb("1e-6")).zeta_nzeros().unique_fmpz()
        right=(h+arb("1e-6")).zeta_nzeros().unique_fmpz()
        assert left is not None and right==left+1
        numbers.append(int(right))
    return numbers


def inside(value,lo,hi):
    if lo < value < hi:
        return True
    assert value < lo or value > hi
    return False


def verify():
    numbers=candidate_numbers()
    assert len(numbers)==69 and numbers[0]==200 and numbers[-2:]==[8297,8359]
    ctx.prec=192
    data=scan(numbers)
    gauge=[]
    curvature=[]
    for n,energy,alpha,curve in data:
        if inside(alpha,arb(".117"),arb(".119")):
            gauge.append(n)
        if inside(curve,arb(".998"),arb("1.002")):
            curvature.append(n)
        print(n,"能标、强耦合、曲率比",energy.str(10),alpha.str(10),curve.str(12))
    assert len(gauge)==23 and curvature==[6514,8035,8359]
    assert set(gauge).intersection(curvature)=={8359}
    prefix={n:1/acb.zeta_zero(n).imag**2 for n in numbers if n<=8192}
    for selected,lo,hi in ((gauge,".3831628",".3831630"),(curvature,".0016815",".0016816")):
        fraction=sum(w for n,w in prefix.items() if n in selected)/sum(prefix.values())
        assert arb(lo)<fraction<arb(hi)
        print("前 8192 枚三模候选中的标准条件概率",fraction.str(12))
    _,_,strong,ratio=data[-1]
    assert arb(".1184112")<strong<arb(".1184213") and arb(".99935101")<ratio<arb(".99935103")
    flavor=joint["flavor"]
    evaluate=flavor["reference_readout"]()
    vub=evaluate(flavor["singlet"](8359)[0])[2]
    reference=evaluate(flavor["singlet"](11661)[0])[2]
    assert arb(".23025")<1-vub/reference<arb(".23026")
    print("首个共同候选的混合读数与参考值",vub.str(12),reference.str(12))


if __name__=="__main__":
    verify()
