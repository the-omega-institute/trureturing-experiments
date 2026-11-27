"""在同一能标处重新约束弱差，在连续能标窗上检验自由荷六截止与中性规范残差。"""
from pathlib import Path
from runpy import run_path

from flint import arb,ctx


def angular_endpoint(gauge,rule,species,p0,p,edges):
    """为一个动量端点 p=M²/M_P² 计算每个截止端点的物种角积分，含 1/(4π)。"""
    assert p>=p0 or p.overlaps(p0)
    if p is p0:
        return [[arb(0) for _ in species] for _ in edges]
    log = (p/p0).log()
    return [[gauge["angular_exposure"](kind,mass,x,p0,p,log,rule)/(4*arb.pi())
             for kind,mass,_ in species] for x in edges]


def hull(parts):
    result = parts[0]
    for part in parts[1:]:
        result = result.union(part)
    return result


def residual_rectangle(previous,current,weights,initial,factors,*,low_end,high_start):
    """从左右动量端点包围两项残差。low_end、high_start 为截止小段的切片界。"""
    rectangle = [[u.union(v) for u,v in zip(left,right,strict=True)]
                 for left,right in zip(previous[:-1],current[1:],strict=True)]
    total = []
    for j in range(2):
        bands = [-sum(b[j]*v for b,v in zip(factors,row,strict=True)) for row in rectangle]
        all_c,low_c,high_c = hull(bands),hull(bands[:low_end]),hull(bands[high_start:])
        total.append(initial[j]+sum(w*(low_c if i==15 else high_c if i==17 else all_c)
                                     for i,w in enumerate(weights)))
    return total


def verify(root):
    def lib(name):
        return run_path(str(root/("ch39-"+name+".py")))

    weights,xmin,xmax,x18,p0,species,initial = lib("fourth-universal-neutral-gauge")["verify"](root)
    lib("neutral-residual-all-cutoffs")["verify"]()
    mp = arb("1.22089e19")
    pmin = (arb("1e13")/mp)**2
    low_fixed = lib("native-rh-ray-pvm")["weak_series"](56)
    def low_uniform(ratio):
        return low_fixed(ratio)+arb(0,(2*ratio*p0/pmin).upper())
    high = lib("full-neutral-free-cutoff")["weak_column"]
    whole,large = lib("outer-shared-distance-kernel")["weak_bands"](low_uniform,high,p0,pmin)
    weak_upper = sum(c*(large if i==15 else whole) for i,c in enumerate(weights)).upper()
    assert weak_upper<arb("-.0152222")
    ctx.prec = 160
    p_orig = (arb("1.176397e14")/mp)**2
    # 相同的截止分割用于所有能标。端点只计算一次，相邻能量小段共享。
    edges = [xmin,4*p0,16*p0,64*p0,256*p0,arb("1e-27"),arb("1e-20"),arb("1e-15"),p_orig,x18,xmax]
    assert all(a<b for a,b in zip(edges,edges[1:]))
    gauge = lib("all-field-signed-crossing")
    rule = gauge["angular_rule"](7,20)
    factors = [(qq-arb(8)*col/3,mix-col) for _,_,(qq,mix,col) in species]
    previous = angular_endpoint(gauge,rule,species,p0,pmin,edges)
    ranges = []
    limits = [(379,682,49,140),(392,709,49,145),(406,740,51,155),(412,785,61,176),
              (409,859,68,210),(417,990,80,256),(428,1269,88,324),(535,1686,103,406)]
    for power in range(27,35):
        energy = arb(10)**(arb(power)/2)
        current = angular_endpoint(gauge,rule,species,p0,(energy/mp)**2,edges)
        total = residual_rectangle(previous,current,weights,initial,factors,low_end=4,high_start=-1)
        p_lo,p_hi,m_lo,m_hi = limits[power-27]
        assert p_lo<total[0]<p_hi and m_lo<total[1]<m_hi
        ranges.append(total)
        print("上端能标的十进指数、两项规范残差的连续外包",power/2,
              [(v.lower().str(12),v.upper().str(12)) for v in total],flush=True)
        previous = current
    print("任意能标的第16层大截止弱差上界",weak_upper.str(13),flush=True)
    assert all(arb("379.304")<p<arb("1685.114") and arb("49.049")<m<arb("405.119") for p,m in ranges)
    return ranges


if __name__=="__main__":
    verify(Path(__file__).resolve().parent)
