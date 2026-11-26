"""计入形状分裂的十六层重 Dirac 后，三支交汇的最小无量纲能标。"""
from pathlib import Path
from runpy import run_path

from flint import acb, arb, ctx


def verify():
    directory = Path(__file__).parent
    shell = run_path(str(directory/"ch39-three-ray-chiral-boundary.py"))
    ctx.prec = 256
    h = acb.zeta_zero(11431).imag
    rows = [[2*x for x in row] for row in shell["model"]["weyl"]]
    qs = shell["model"]["charges"]
    levels = {q:sorted((v.real for v in shell["laplacian"](q,h).eig()),key=float)[3:] for q in qs}
    a = [3*(1-arb(".23129"))*arb("127.930")/5, arb(".23129")*arb("127.930"), 1/arb(".1180")]
    b1,b2,b3 = arb(41)/10,-arb(19)/6,-7
    # Certified by ch39-full-shell-haar-grid.py for the actual two-step face.
    r = arb("60.29235").union(arb("60.29245"))/60
    s = arb(".9987274").union(arb(".9987286"))
    normal = ((b3-s*b2)/(r*b2-b1), (s*b1-r*b3)/(r*b2-b1), arb(1))
    target = sum(x*y for x,y in zip(normal,a))
    projection = {q:sum(x*y for x,y in zip(normal,row)) for q,row in zip(qs,rows)}
    doublets = (1,3)
    assert all((c < 0 if q in doublets else c > 0) for q,c in projection.items())
    sums = {q:sum(x.log() for x in levels[q])/2 for q in qs}
    log_onset = ((2*arb.pi()*target+sum(projection[q]*sums[q] for q in doublets)) /
                 (16*sum(projection[q] for q in doublets)))
    onset = log_onset.exp()
    assert onset**2 > max(x for q in qs for x in levels[q])
    theta = shell["geo"]["theta"]
    native = theta(804)/theta(12)
    assert arb("9.659") < onset < arb("9.661") and onset/native > arb("1.1317")
    native_gap = sum(projection[q]*(16*native.log()-sums[q]) for q in doublets)/(2*arb.pi())-target
    assert native**2 > max(x for q in doublets for x in levels[q]) and native_gap > arb(".775")

    def run_match(tower, ratio):
        at = [(a[i]-sum(row[i]*tower[q] for q,row in zip(qs,rows))/(2*arb.pi())) for i in range(3)]
        clock = (at[0]-r*at[1])/(b1-r*b2)
        strength = at[1]-b2*clock
        scale = arb("91.1876")*(2*arb.pi()*clock).exp()
        assert strength > 0 and scale/ratio > arb("1e10")
        assert (at[2]-s*strength-b3*clock).contains(0)
        return strength,scale

    exposure = {q:16*log_onset-sums[q] if q in doublets else arb(0) for q in qs}
    min_lift = onset**2-min(levels[q][0] for q in qs if q not in doublets)
    assert arb("83.62") < min_lift < arb("83.66")
    k,scale = run_match(exposure,onset)
    print("首次容许能标比、几何能标的投影缺口、临界单态共同平方质量提升",onset,native_gap,min_lift)
    print("临界处的共同归一、交汇能标 GeV",k,scale)
    ratio = arb(10)
    def split(lift):
        return {q:sum((ratio/(v+(0 if q in doublets else lift)).sqrt()).log().max(0)
                      for v in levels[q]) for q in qs}
    def residual(lift):
        return sum(projection[q]*value for q,value in split(lift).items())/(2*arb.pi())-target
    lo,hi = arb("77.54"),arb("77.58")
    assert residual(lo)>0 and residual(hi)<0
    assert any(levels[q][0]+hi < ratio**2 for q in qs if q not in doublets)
    k,scale = run_match(split(lo.union(hi)),ratio)
    assert arb("5.96e14") < scale < arb("6.01e14")
    print("能标比 10 的唯一共同单态提升区间、共同归一与能标",lo.union(hi),k,scale)


if __name__ == "__main__":
    verify()
