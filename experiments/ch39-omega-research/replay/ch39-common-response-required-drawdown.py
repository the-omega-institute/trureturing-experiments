"""两枚实际编号在任意复味空间上匹配顶底质量与两项引力所需的能序回落。"""
from pathlib import Path
from runpy import run_path

from flint import acb, arb, ctx


def field_bounds(rows):
    pi = arb.pi()
    k = -sum(a for _,a,_ in rows)/(2*pi)
    w = -sum(g for g,_,_ in rows)/(8*pi)
    assert k>0 and w>0 and all(z<1 for _,_,z in rows)
    massive = [(g,a,z) for g,a,z in rows if not z.is_zero()]
    assert all(z>0 for _,_,z in massive)
    err_k = sum(abs(a)*z*(2-z.log()) for _,a,z in massive)/(2*pi)
    err_w = sum(abs(g)*z for g,_,z in massive)/(4*pi)
    return k,w,err_k,err_w


def check(upper_drop,alpha,weights,rows):
    # 只用同表顶的下端和底的上端；其余味观测不进入验算。
    lower_top,upper_bottom = arb(".8573"),arb(".01410")
    d = upper_drop
    x_min = (lower_top-d)**2
    small = (upper_bottom/alpha+d)**2
    assert lower_top>d and x_min>0
    # 记 X=x19，则每个荷六位置 x_j≤(√X+d)²，max(1,x_j)≤scale·X。
    scale = max((1/x_min).upper(),((1+d/(lower_top-d))**2).upper())
    total = sum(abs(c) for c in weights)
    low_total = sum(abs(c) for c in weights[:16])
    p,b,c = weights[16],-weights[17],weights[18]
    assert b>p>0 and c>0
    k,w,err_k,err_w = field_bounds(rows)
    ek = total*err_k*scale/k
    ew = total*err_w*scale/w
    # K 方程给 bZ−pY>0，并给其除以 X 后的上界。
    assert c>ek and 1/k>low_total*small
    linear = c+ek+(1/k+low_total*small)/x_min
    # W 方程给 (bZ²−pY²)/X² 的下界。
    square = c-ew/x_min-(low_total*small**2+arb("1e-122")/w)/x_min**2
    assert square>0
    gap = ((b-p)*square).sqrt()-linear
    assert gap>0
    return gap


def verify(root):
    native = run_path(str(root/"ch39-free-four-site-native-positive-spectrum.py"))
    shell = run_path(str(root/"ch39-three-ray-chiral-boundary.py"))
    sorted_eigen = run_path(str(root/"ch39-bounded-common-rational-amplitude.py"))["ordered"]
    pivots = run_path(str(root/"ch39-native-rh-ray-pvm.py"))
    gravity = run_path(str(root/"ch39-three-index-free-atoms.py"))
    numbers = (2_306_888_100_085,2_306_888_100_086)
    _,_,fields = gravity["broken_fields"](root,numbers)
    ctx.prec = 320
    for n,rows,drop in zip(numbers,fields,(".013",".050"),strict=True):
        _,down = native["native"](root,n)
        tail = down[-2]+down[-1]
        trace = tail.trace().real
        # 两个秩一块之和至多秩二，非零本征值由迹和迹平方恢复。
        gap_squared = 2*(tail*tail).trace().real-trace**2
        assert gap_squared>0
        alpha = (trace-gap_squared.sqrt())/2
        assert alpha>arb(".90")
        hh = acb.zeta_zero(n).imag
        up_levels,down_levels = [sorted_eigen(shell,q,hh) for q in (-4,2)]
        _,six,projections,_ = next(t for t in pivots["shell_weights"](root,n) if t[0]==6)
        assert six[15] < down_levels[17] and up_levels[-1] < six[-1]
        weights = [v.real for v in projections[0]]
        delta = arb(drop)
        gap = check(delta,alpha,weights,rows)
        print("编号、下型末两投影次大本征值、必需回落的严格下界、无量纲矛盾余量",
              n,alpha.str(12),delta.str(5),gap.str(11),flush=True)


if __name__=="__main__":
    verify(Path(__file__).resolve().parent)
