"""按全谱逆高度平方概率，对固定真空的 105 圈十二面因子给出统一均方下界。"""
from pathlib import Path
from runpy import run_path
import sys

from flint import acb, arb, arb_mat, ctx


# 放宽到只保留归一性及两个有理高阶下界。所列阶在最优分布中严格超出下界。
FREE = (1,3,6,8,11,12,13,17,19,20,22,23,24,25,26,27,28,29,31,33,34,35,37,40,
        42,43,44,46,47,50,59,62,64,65,66,68,70,71,75,76,77,79,80,81,82,83,85,86,88,89,96,97)


def verify(root):
    source = run_path(str(root/'ch39-three-mode-arithmetic-mass.py'))
    native = run_path(str(root/'ch39-mincut-real-triangle-tomography.py'))
    fixed = run_path(str(root/'ch39-height-independent-flavor-vacuum.py'))
    lib = run_path(str(root/'ch39-mincut-global-winding.py'))
    angles = run_path(str(root/'ch39-mincut-reference-ckm-action.py'))
    flavor = run_path(str(root/'ch39-chiral-neutral-flavor-phase.py'))
    ctx.prec, ctx.cap = 224, 3
    total = source['time_data']()[0]
    assert arb('.023')<total<arb('.024')
    shell = source['shell']
    x=arb('-.01046')
    y=x-angles['reference_angles'](flavor)[0]
    tri=lib['terms'](shell)
    fixed['required_high_weights'](lib,tri,x,y)
    coeff=[lib['coeff'](tri,q,arb(0))[2].real/12 for q in range(1,106)]
    u,v = [native['unit'](label) for label in ((13872,23),(-26877,115))]
    assert all((c-((u**q).real+(v**q).real)/2).contains(0) for q,c in enumerate(coeff,1))
    geo=shell['geo']
    frequency=geo['b']-geo['a']
    zeros=acb.zeta_zeros(1,512)
    assert len(zeros)==512 and all(z.real==arb(1)/2 for z in zeros)
    assert all(a.imag<b.imag for a,b in zip(zeros,zeros[1:]))
    heights=[z.imag for z in zeros]
    raw_mass=sum(1/h**2 for h in heights)
    design=arb_mat([[c*(q*frequency*h).cos()/h for q,c in enumerate(coeff,1)] for h in heights])
    gram=design.transpose()*design/total
    # 线性目标和互补条件均在全谱归一的 512 项半正定二次型上核验。
    floor={104:arb('.006589'),105:arb('.068674')}
    free=set(FREE)
    assert not free.intersection(floor)
    kkt=arb_mat([[gram[i-1,j-1] for j in FREE]+[-1] for i in FREE]
                +[[1]*len(FREE)+[0]])
    rhs=arb_mat([[-sum(t*gram[i-1,j-1] for j,t in floor.items())] for i in FREE]
                +[[1-sum(floor.values())]])
    solution=kkt.solve(rhs).entries()
    multiplier=solution[-1]
    weights={**floor,**dict(zip(FREE,solution[:-1],strict=True))}
    assert min(solution[:-1])>0 and sum(weights.values()).contains(1)
    residual=gram*arb_mat([[weights.get(q,0)] for q in range(1,106)])
    for q in range(1,106):
        if q in free:
            assert (residual[q-1,0]-multiplier).contains(0)
        else:
            assert residual[q-1,0]>multiplier
    minimum=sum(w*residual[q-1,0] for q,w in weights.items())
    assert minimum.rad()<arb('1e-35')
    print('全谱总权、首 512 枚归一质量、末枚高度',total.str(18),(raw_mass/total).str(15),heights[-1].str(16))
    print('精确放宽最小值、自由支撑数、均方根',minimum.str(20),len(FREE),minimum.sqrt().str(15))
    assert minimum>arb('3.10017039624786e-5')
    assert minimum.sqrt()>arb('.005567917381')


if __name__=='__main__':
    verify(Path(sys.argv[1]) if len(sys.argv)>1 else Path(__file__).resolve().parent)
