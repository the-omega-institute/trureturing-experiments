"""求右弱环重复算术相位时，选择 11431 的最少端口次数。"""
from flint import acb, arb, ctx


def angle(index):
    return 2*(1/(2*acb.zeta_zero(index).imag)).atan()


def count(height):
    n = height.zeta_nzeros().unique_fmpz()
    assert n is not None
    return int(n)


def near_resonances(power, width):
    first_angle = angle(1)
    spacing = 2*arb.pi()/3
    candidates = set()
    intervals = 0
    m = 1
    while m*spacing-width <= power*first_angle:
        if m%3:
            center = m*spacing
            lo = 1/(2*((center+width)/(2*power)).tan())
            hi = 1/(2*((center-width)/(2*power)).tan())
            for n in range(count(lo.lower())+1,count(hi.upper())+1):
                candidates.add((n,m))
            intervals += 1
        m += 1
    return intervals, candidates


def verify():
    ctx.prec = 128
    target = 11431
    t0,theta,t2 = [angle(n) for n in (target-1,target,target+1)]
    resonance = 2*arb.pi()/3
    lower = 2*resonance/(t0+theta)
    upper = 2*resonance/(theta+t2)
    assert 23225 < lower < 23226 < upper < 23227
    power = 23226
    delta = abs(power*theta-resonance)
    width = arb("0.00016")
    assert delta < width < resonance/2
    bands, pairs = near_resonances(power,width)
    readings = sorted([(abs(power*angle(n)-m*resonance),n,m) for n,m in pairs],
                      key=lambda row:float(row[0].mid()))
    assert bands == 523 and [(n,m) for _,n,m in readings]==[(11431,1),(1144,7),(11432,1)]
    assert all(left[0] < right[0] for left,right in zip(readings,readings[1:]))
    best,next_best = readings[0][0],readings[1][0]
    assert next_best < width and 2*power*arb("7.5e-10") < next_best-best
    gap = 2*(best.cos()-next_best.cos())
    assert gap > arb("3.04e-9")
    print("first integer interval",lower.str(13),upper.str(13),"power",power,"bands",bands,
          "contenders", [(n,m,d.str(8)) for d,n,m in readings[:8]],flush=True)
    print("path eigenvalue gap / g",gap.str(12),flush=True)


if __name__=="__main__":
    verify()
