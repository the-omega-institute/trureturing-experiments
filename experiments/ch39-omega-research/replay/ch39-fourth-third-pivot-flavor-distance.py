"""对所有第三枢轴五式解，界定其与第四壳已证完整味见证的射影算子距离。"""
from pathlib import Path
from runpy import run_path

from flint import arb,ctx


def angles(flavor,source):
    low,top = source[0].sqrt(),source[-1].sqrt()
    initial = [((x.sqrt()-low)/(top-low)).sqrt().asin() for x in flavor]
    return low,top,initial


def row_bounds(flavor,source,phi):
    low,top,initial = angles(flavor,source)
    bounds = []
    for s in initial:
        lo = max(arb(0),(s-phi).lower())
        hi = min((arb.pi()/2).upper(),(s+phi).upper())
        bottom = low+(top-low)*lo.sin()**2
        upper = top if (s+phi).upper()>=(arb.pi()/2).lower() else low+(top-low)*hi.sin()**2
        bounds.append((bottom.lower(),upper.upper()))
    return bounds


def thermal_band(bounds,weights,slopes,errors,mix):
    k,w = slopes
    total = arb(0)
    for (l,r),c in zip(bounds,weights,strict=True):
        def value(f2):
            return f2*(k+mix*w*f2)
        left,right = l**2,r**2
        val = value(left).union(value(right))
        if mix>0:
            crit = -k/(2*mix*w)
            if not (crit<left or right<crit):
                val = val.union(value(crit))
        total += c*val
    return total+arb(0,(errors[0]+abs(mix)*errors[1]).upper())


def verify(root):
    def lib(name):
        return run_path(str(root/("ch39-"+name+".py")))
    n = 2_306_888_100_085
    xflavor,_ = lib("fourth-full-ckm-gravity")["verify"](root,return_cutoffs=True)
    _,source = lib("adjacent-ordered-charge-six-cutoffs")["verify"](root,refine_node=n+1)[2]
    native = lib("native-rh-ray-pvm")
    _,_,pivots,_ = next(row for row in native["shell_weights"](root,n+4) if row[0]==6)
    weights = [q.real for q in pivots[2]]
    moments,heat,(fields,) = lib("three-index-free-atoms")["broken_fields"](root,[n+4])
    ctx.prec = 192
    pi = arb.pi()
    slopes = [sum(a for _,a,_ in fields)/(2*pi),-sum(g for g,_,_ in fields)/(8*pi)]
    xmax = source[-1]
    tv = sum(map(abs,weights))
    errors = [tv*xmax*sum(abs(a)*(1-heat(2,z/xmax)) for _,a,z in fields)/(2*pi),
              tv*xmax**2*sum(abs(g)*(arb(1)/2-heat(3,z/xmax)) for g,_,z in fields)/(4*pi)]
    assert slopes[0]<0<slopes[1]
    # 第16层虽在电弱尺度，完整热迹列在这里同样用原始破缺场核。
    cols = [moments(fields,heat,x) for x in xflavor]
    original = [sum(c*v[j] for c,v in zip(weights,cols,strict=True)) for j in range(2)]
    print("原味见证在第三枢轴的 K、W",*[v.str(14) for v in original],flush=True)
    assert arb("172.07380")<original[0]<arb("172.07381")
    assert arb("-119.70868")<original[1]<arb("-119.70866")
    phi = arb(".36")
    mix = arb(19)
    reachable = thermal_band(row_bounds(xflavor,source,phi),weights,slopes,errors,mix)
    target = 1+mix*arb("1e-122")
    assert reachable.upper()<arb("-1.38668")<target
    cubic = 2*(phi/2).sin()
    generator = 2*(phi/6).sin()
    assert cubic>arb(".3580") and generator>arb(".11992")
    print("0.36 弧度内 K+19W 上界、所有解的立方与生成元距离严格下界",
          *[v.str(14) for v in (reachable.upper(),cubic,generator)],flush=True)

    witness = lib("fourth-third-pivot-five")["verify"](root)
    low,high,old = angles(xflavor,source)
    pos = ((witness[18]-low)/(high-low)).sqrt().asin()
    gap = (old[18]-pos).lower()
    assert gap>arb("1.0991")
    box_cubic = 2*(gap/2).sin()
    box_gen = 2*(gap/6).sin()
    assert box_cubic>arb("1.0446") and box_gen>arb(".3643")
    print("五式解盒第19层角差、立方及所有立方根的距离下界",
          *[v.str(14) for v in (gap,box_cubic,box_gen)],flush=True)
    return original,reachable,gap


if __name__=="__main__":
    verify(Path(__file__).resolve().parent)
