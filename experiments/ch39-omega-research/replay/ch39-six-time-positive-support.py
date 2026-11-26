"""六时间全支撑逼近及两枚实际零点不能同时消去的定量证书。"""
from pathlib import Path
from runpy import run_path
import sys

from flint import acb, arb, arb_mat, ctx


PAIR = (142_270_715_263, 143_921_351_258_535)
SCALE = 10**11
MIX = (74806494524, 25193505476)
GRADIENT = (-223155029, -36853801, -67588124, -28770032)
RIVALS = ((1389123873, 14439811883, 3195483964),
          (1696866033, 9230799675, 782619241),
          (-1691586351, -27495589798, 66923703),
          (1566585886, 4910246391, 122609047))


def time_column(height, scale):
    u = (scale/height)**2
    return [arb(1), *(1/(1+j*j*u) for j in range(1,7))]


def verify(root):
    def read(name):
        return run_path(str(root/('ch39-'+name+'.py')))
    exact = read('six-time-proton-exact-zero')
    support = exact['SUPPORT']
    assert not set(PAIR).intersection(support)
    heights, probabilities = read('six-time-proton-small-mean')['six_time_support'](root, support)
    ctx.prec = 288
    scale = read('six-time-height-bands')['scale']
    matrix = arb_mat([time_column(heights[n],scale) for n in support]).transpose()
    inverse = matrix.inv()
    budgets = [p/sum(abs(inverse[i,j]) for j in range(7))
               for i,p in enumerate(probabilities)]
    budget = min(z.lower() for z in budgets)
    assert budget>arb('2.43494e-15')
    # 对任意额外字典及其先验，七个时间坐标都位于 [0,1]。
    print('任意额外先验的充分掺入量上界（各点约束的最小值）',budget.str(17))

    shell = read('three-ray-chiral-boundary')
    lib, fixed = read('mincut-global-winding'), read('height-independent-flavor-vacuum')
    flavor, angles = read('chiral-neutral-flavor-phase'),read('mincut-reference-ckm-action')
    physical = read('first-three-mode-observables')
    model, initial = physical['model_inputs']()
    chiral = physical['joint']['chiral']
    ctx.prec = 192
    gaps=[]
    for n in PAIR:
        height = acb.zeta_zero(n).imag
        clocks = shell['geo']['phases'](height)
        gaps.append((clocks[724]/clocks[804]).arg())
        rows = chiral['source']['spectral_supplier'](model)(height)
        layers = {q:chiral['ordered']([v for v,z in rows if z==q]) for q in (0,*model['charges'])}
        _,alpha,curve = physical['layer_observables'](model,initial,layers)
        assert arb('.117')<alpha<arb('.119') and arb('.998')<curve<arb('1.002')
        print('额外编号、钟角差、强耦合和曲率比',n,*(v.str(15) for v in (gaps[-1],alpha,curve)))
    ctx.prec=224
    x=arb('-.01046')
    y=x-angles['reference_angles'](flavor)[0]
    tri=lib['terms'](shell)
    mix=[arb(z)/SCALE for z in MIX]
    gradients=[arb(z)/SCALE for z in GRADIENT]
    rivals=[[arb(v)/SCALE for v in row] for row in RIVALS]
    assert sum(MIX)==SCALE and all(r>0 for _,_,r in rivals)
    margins=[]
    for q in range(1,106):
        terms=lib['coeff'](tri,q,arb(0))
        reading=sum(v*(terms[2]*acb(0,q*d).exp()).real/12 for v,d in zip(mix,gaps,strict=True))
        base=lib['energy'](q,terms,x,y)
        compares=sum(r*(lib['energy'](q,terms,a,b)-base) for a,b,r in rivals)
        gradient=sum(h*g for h,g in zip(gradients,fixed['gradient_columns'](lib,tri,q,x,y),strict=True))
        margins.append(reading-gradient-compares)
    bound=min(z.lower() for z in margins)
    constant=bound**2/sum(v**2 for v in mix)
    assert bound>arb('.00456315') and constant>arb('.0000334188')
    print('两点凸组合的逐阶严格下界、两概率共同下界的均方系数',bound.str(15),constant.str(15))


if __name__=='__main__':
    verify(Path(sys.argv[1]) if len(sys.argv)>1 else Path(__file__).resolve().parent)
