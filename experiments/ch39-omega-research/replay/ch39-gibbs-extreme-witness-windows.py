"""两个严格三窗中的作用间隔与必要高度段的逆平方总量。"""
from pathlib import Path
from runpy import run_path
import sys

from flint import acb, arb, ctx


root = Path(sys.argv[1]) if len(sys.argv)>1 else Path(__file__).resolve().parent
old = run_path(str(root/"ch39-triangle-gibbs-six-times.py"))
windows = run_path(str(root/"ch39-three-window-extreme-cp.py"))


def witnesses():
    coords = ((".07174503",".07143683"),(".16070319","-.05199108"))
    windows["open_windows"](coords)
    ctx.prec=192
    counts, coeff, base = old["harmonic_data"]()
    for i,(x,y) in enumerate(coords):
        r=arb(2)**-42
        u,v=[acb(0,arb(t)+arb(0,r)).exp() for t in (x,y)]
        value=sum(n-(coeff[a,b,j]*u**a*v**b).real for (a,b,j),n in counts.items())-base
        # 真实零点高度 >10000 时的软边误差；每个软化三角至多贡献 2δ。
        err=sum(2*j*n for (_,_,j),n in counts.items())*arb(23)/25*arb(2)**-20/10000**2
        value+=arb(0,err)
        if i==0:
            assert value<arb(".15313021")
        else:
            assert value>arb(".42630569")
        print("开窗中的作用增量",value.str(15))
    lo,hi=arb(41339000000),arb(41399000000)
    def main(t):
        return t/(2*arb.pi())*((t/(2*arb.pi())).log()-1)+arb(7)/8
    N=main(hi)-main(lo)+arb(".28")*(hi.log()+lo.log())
    prior=N/lo**2
    assert prior<arb("1.263317e-13")
    low_cap=arb("4.676839e-11")
    for number, floor, edge in ((117930158,"98714.51",".1531305"),(71363114,"12158.64",".426303")):
        h, energy=old["record"](number,(counts,coeff,base))
        local=(arb(".10448811")/(low_cap*prior*h**2)).log()/abs(energy-arb(edge))
        assert local > arb(floor)
        print("窄段赔率必要逆温度绝对值",number,local.str(13))
    print("必需窄段未归一化逆平方总量上界",prior.str(15))


if __name__=="__main__":
    witnesses()
