"""共同正电导二步算子的荷间重层质量和上界，以及绝对质量归一化的简并。"""
from fractions import Fraction
from pathlib import Path
from runpy import run_path
import sys

from flint import acb, arb, ctx, fmpq_mat


NUMBER=2_306_888_100_085


def verify(root):
    shell=run_path(str(root/"ch39-three-ray-chiral-boundary.py"))
    graph=shell["geo"]
    edges=[(i,j) for i,j,_,_ in graph["edges"]]
    assert len(edges)==135 and len(graph["sites"])==19
    # H>0；以总电导为1，在任意共享电导上精确极小化不依赖荷的二阶谱矩。
    incidence=fmpq_mat(19,135)
    for k,(i,j) in enumerate(edges):
        incidence[i,k]=incidence[j,k]=1
    h=incidence.transpose()*incidence
    for i in range(135):
        h[i,i]+=2
    ones=fmpq_mat([[1] for _ in edges])
    v=h.solve(ones)
    z=sum(v.entries())
    assert min(v.entries())>0
    sigma=81/z-18
    ctx.prec=256
    assert int(sigma.p)==57528 and int(sigma.q)==144905
    s=arb(int(sigma.p))/int(sigma.q)
    C=(s+18).sqrt()*(s+1).sqrt()/s
    upper=arb(19)/16*C
    print("共用135条电导的最小归一化有符号迹",sigma,s.str(18))
    print("任意荷两完整质量和、十六重层质量和之比的上界",C.str(18),upper.str(18))
    height=acb.zeta_zero(NUMBER).imag
    delta=arb(23)/25*arb(2)**-20/(height**2+shell["soft"]["point"]["study"]["offset"])
    c=arb(15)-delta/3
    bands={}
    for q in (1,2):
        lap=shell["laplacian"](q,height)
        eig=lap.eig(multiple=True)
        assert len(eig)==19 and all(z.imag.contains(0) and z.real>0 for z in eig)
        masses=[z.real*abs(z.real-c) for z in eig]
        ordered=sorted(masses,key=float)
        assert all(a<b for a,b in zip(ordered,ordered[1:]))
        bands[q]=sum(ordered[3:])
    native=bands[2]/bands[1]
    physical=32*native
    print("真实零点荷一、荷二的十六层和、原比与乘32后的比",
          bands[1].str(18),bands[2].str(18),native.str(18),physical.str(18))
    print("保持两塔原形时，共同正电导可生成的任何偶倍率上界",(upper/native).str(18))
    assert arb("30.6226")<physical<arb("30.6227") and upper<arb("15.163877")
    assert upper/native<arb("15.846") and physical/upper>arb("2.019")
    print("即使边权自由仍需另加的偶奇系数比下界",(physical/upper).str(18))
    # t 倍的电导给 t² 倍的二步质量。全部原模型质量在 d->td 后逐层相同。
    d=Fraction(724)
    t=Fraction(2)
    assert t*t/(t*d)**2==1/d**2
    print("电导单位、分母和电导平方/分母平方的两组完全相同质量",1,d,t,t*d,1/d**2)


if __name__=="__main__":
    verify(Path(sys.argv[1]) if len(sys.argv)>1 else Path(__file__).resolve().parent)
