"""同一归一截止同时匹配牛顿及独立弱差所需的带宽与源重叠。"""
from pathlib import Path
from runpy import run_path
from flint import arb,ctx


def verify(root):
    ctx.prec=192
    _,fields=run_path(str(root/"ch39-three-index-free-atoms.py"))["spectra"](root)
    largest=sum(abs(a) for _,a,_ in fields[0])
    assert all((sum(abs(a) for _,a,_ in group)-largest).contains(0) for group in fields)
    # 分正负号给更紧的 K 上界。积分核 E2 非负且至多一。
    plus = sum(a for _,a,_ in fields[0] if a>0)
    minus = -sum(a for _,a,_ in fields[0] if a<0)
    for group in fields:
        assert all(a.is_zero() or a>0 or a<0 for _,a,_ in group)
        assert (sum(a for _,a,_ in group if a>0)-plus).contains(0)
        assert (-sum(a for _,a,_ in group if a<0)-minus).contains(0)
    b=arb("3844.11")
    planck=arb("1.22089e19")
    high=planck*(2*arb.pi()/((1+b)*plus+b*minus)).sqrt()
    largest_weak=arb("2662.44")
    ratio=high/largest_weak
    necessary=(2*arb.pi()*arb("9.9443061e33")/max(plus,minus)).sqrt()
    print("热核牛顿正负系数、最大截止下界 GeV、所需能量比、总变差乘带宽的下界",
          plus.str(10),minus.str(10),high.lower().str(14),ratio.lower().str(13),
          necessary.lower().str(13))
    native=run_path(str(root/"ch39-native-rh-ray-pvm.py"))
    for q,levels,_,_ in native["shell_weights"](root):
        bandwidth=(levels[-1]/levels[0]).sqrt()
        tv=necessary/bandwidth
        print("荷、原谱平方质量读法的带宽与总变差下界",q,
              bandwidth.upper().str(11),tv.lower().str(12))
        assert bandwidth<ratio and tv>arb("1e12")
    for group in fields:
        heavy=[z for _,_,z in group[4:]]
        bandwidth=(max(z.upper() for z in heavy)/min(z.lower() for z in heavy)).sqrt()
        print("该零点重谱的能量带宽上界",bandwidth.upper().str(12))
        assert bandwidth<ratio
    assert high>arb("3.4741e16") and ratio>arb("1.3048e13")
    assert necessary>arb("2.3887e16")


if __name__=="__main__":
    verify(Path(__file__).resolve().parent)
