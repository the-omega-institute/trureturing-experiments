"""以真实源谱的最低截止和弱差核两端重算第八壳色阈值的严格上界。"""
from inspect import getclosurevars
from pathlib import Path
from runpy import run_path

from flint import arb, ctx


def weak_limits(lib, root):
    _, ordered = lib("adjacent-ordered-charge-six-cutoffs")["verify"](
        root, refine_node=2_306_888_100_086)[2]
    ctx.prec = 256
    p0 = (arb("91.1876")/arb("1.22089e19"))**2
    p1 = (arb("1.176397e14")/arb("1.22089e19"))**2
    inv_max = p0/ordered[0]
    assert inv_max < arb(".784")
    series = getclosurevars(lib("native-rh-ray-pvm")["weak_series"](56)).nonlocals
    coeff = series["coefficients"]
    # P'(s)>0，s=p0/x；一次取八个闭段覆盖包括大截止极限 s=0 的全区间。
    slopes = []
    for j in range(8):
        s = arb(j).union(arb(j+1))*arb(".784")/8
        value = arb(0)
        for n in range(len(coeff),0,-1):
            value = value*s + n*coeff[n-1]
        slopes.append(value.lower())
    assert min(slopes) > arb(".01567")
    value = arb(0)
    for a in reversed(coeff):
        value = value*inv_max+a
    max_correction = inv_max*value
    join = arb(10)**10
    radius = (series["remainder"] + 2*join*series["p_ratio"]).upper()
    moderate_lo = -100*(series["reference"] + max_correction + radius)
    moderate_hi = -100*(series["reference"] - radius)
    # weak_column 的余量随 x 递减，故连接点的界覆盖它右侧的全部截止。
    large = lib("full-neutral-free-cutoff")["weak_column"](join*p0,p0,p1)
    lower, upper = arb("5.37726334919"),arb("7.425657721")
    assert lower < moderate_lo and moderate_hi < upper and lower < -100*large < upper
    print("源谱最小对数截止、弱差级数导数下界、负百倍弱差范围",
          (-inv_max.log()).str(16), min(slopes).str(13), lower.str(13), upper.str(13), flush=True)
    return lower, upper


def verify(root):
    def lib(name):
        return run_path(str(root/("ch39-"+name+".py")))

    lo,hi = weak_limits(lib, root)
    support = lib("cross-shell-pivot-support")
    core = lib("fourth-third-pivot-global-kinetic-obstruction")
    arc = lib("eighth-shortest-arc-pivot-support")
    rational = lib("eighth-rational-arc-signs")
    raw = support["projected_tables"](core, root, arc["NUMBER"], None)
    scale = core["SCALE"]
    lower = int((lo*scale).lower().floor().unique_fmpz())
    upper = int((hi*scale).upper().ceil().unique_fmpz())
    tables = []
    for (offset,bands,unit), d in zip(raw,support["DIRECTIONS"],strict=True):
        if d == (0,0,-1,0):
            assert unit == scale
            bands = [(max(l,lower),min(r,upper)) for l,r in bands]
            assert all(l <= r for l,r in bands)
        tables.append((offset - d[3]*arb(2)*scale/10,bands,unit))
    rows = arc["polynomials"](root,include_single_edge=False)[0][1]
    left,right = ".73613582",".73613585"
    for t,sgn in ((left,-1),(right,1)):
        w = rational["rational_weights"](rows,arb(t))
        a,b = arc["score"](tables,w)
        assert a>0 if sgn>0 else b<0
        print("参数与修正支撑指标的整数外包（单位千万分之一）", t, int(a), int(b),flush=True)
    rational["crossing"](arc,tables,rows,left,right,(0,7,16),1)
    rational["check_sign"](arc,tables,rows,right,".75",1)


if __name__ == "__main__":
    verify(Path(__file__).resolve().parent)
