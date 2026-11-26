"""重建十二面的整数标签，排除全角恒零并给出 105 圈模型的均方下界。"""
from collections import Counter
from fractions import Fraction
from pathlib import Path
from runpy import run_path
import sys

from flint import arb, ctx


def squared_phase_trace(label, shell_radius):
    a, b = label
    return Fraction(2*(a*a-shell_radius*b*b), a*a+shell_radius*b*b)


def verify(root):
    native = run_path(str(root/'ch39-mincut-real-triangle-tomography.py'))
    fixed = run_path(str(root/'ch39-height-independent-flavor-vacuum.py'))
    winding = run_path(str(root/'ch39-mincut-global-winding.py'))
    shell = run_path(str(root/'ch39-three-ray-chiral-boundary.py'))
    angles = run_path(str(root/'ch39-mincut-reference-ckm-action.py'))
    flavor = run_path(str(root/'ch39-chiral-neutral-flavor-phase.py'))
    ctx.prec = 224
    counts = Counter(label for k,label in native['enumerate_faces']().values() if k == (1,-1))
    left, right = (13872,23), (-26877,115)
    conjugate = native['conjugate']
    assert counts == Counter({label:3 for label in (left,conjugate(left),right,conjugate(right))})
    ratio_labels = [native['product'](left, label) for label in (right, conjugate(right))]
    traces = [squared_phase_trace(label,native['R']) for label in ratio_labels]
    assert traces == [Fraction(1899063107,960768307), Fraction(64648907,34388347)]
    assert all(z.denominator>1 for z in traces)
    print('乘积与商的平方之有理迹', *traces)

    u, v = [native['unit'](label) for label in (left,right)]
    C = {q:((u**q).real+(v**q).real)/2 for q in range(1,106)}
    tri = winding['terms'](shell)
    assert len(tri[1,-1])==12 and all(soft==0 for _,soft in tri[1,-1])
    for q, coeff in C.items():
        original = winding['coeff'](tri,q,arb(0))[2]/12
        assert not coeff.contains(0) and (original-coeff).contains(0)
    baseline = 1/(2*sum(1/z**2 for z in C.values()))
    assert arb('.00000198398919756689')<baseline<arb('.00000198398919756690')
    x=arb('-.01046')
    y=x-angles['reference_angles'](flavor)[0]
    fixed['required_high_weights'](winding,tri,x,y)
    floors={104:arb('.006589'),105:arb('.068674')}
    bulk=sum(1/C[q]**2 for q in range(1,104))
    remainder=1-sum(floors.values())
    # 固定高阶份额后，其余各阶按 C_q^{-2} 分配，使均方最小。
    # 在两个高阶下界处两条偏导均正，凸性给出整个可行象限的下界。
    assert all(C[q]**2*f>remainder/bulk for q,f in floors.items())
    constrained=(remainder**2/bulk+sum((C[q]*f)**2 for q,f in floors.items()))/2
    assert arb('.00011962850134974')<constrained<arb('.00011962850134976')
    weights = fixed['minimum_order_weights'](winding,tri,x,y)
    example=sum((w*C[q])**2 for q,w in weights.items())/2
    assert example>constrained
    print('仅归一时的精确最小均方、纳入强制高阶份额的均方下界、固定真空示例的均方',
          *(z.str(17) for z in (baseline,constrained,example)))
    print('对应强制下界的均方根', constrained.sqrt().str(15))


if __name__=='__main__':
    verify(Path(sys.argv[1]) if len(sys.argv)>1 else Path(__file__).resolve().parent)
