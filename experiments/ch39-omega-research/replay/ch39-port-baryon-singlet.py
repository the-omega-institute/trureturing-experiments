"""全谱选点与同一两臂电流的轮换保护之间的定量阻碍。"""
from pathlib import Path
from runpy import run_path

from flint import arb, acb, ctx

root = Path(__file__).with_name("ch39-pure-port-direction-window.py")
direction = run_path(str(root))
point = run_path(str(root.with_name("ch39-port-coupling-competition.py")))


def verify():
    ctx.prec = 192
    constraints = direction["inequalities"]()
    stack = [(arb("2.220"),arb("2.242"),arb(".892"),arb(".920"),0)]
    count = 0
    while stack:
        lo,hi,left,right,depth = stack.pop()
        alpha,phi = lo.union(hi),left.union(right)
        x = alpha.cos()
        p = (1+x)/2
        cosine = (5*phi).cos()
        total = p+(1-p)*cosine
        difference = total/(1-cosine)
        count += 1
        if total > arb(".005") and difference > arb(".004"):
            continue
        y,z = alpha.sin()*phi.cos(),-alpha.sin()*phi.sin()
        if any(a+b*x+c*y+e*z < 0 for a,b,c,e in constraints):
            continue
        assert depth<22
        if hi-lo>right-left:
            mid=(hi+lo)/2
            stack.extend([(lo,mid,left,right,depth+1),(mid,hi,left,right,depth+1)])
        else:
            mid=(left+right)/2
            stack.extend([(lo,hi,left,mid,depth+1),(lo,hi,mid,right,depth+1)])
    vec = point["prepare"]()
    phase = (vec[1,0]/abs(vec[1,0]))**5
    old = arb(5)/26
    altered = -phase.real/(1-phase.real)
    transfer = old-altered
    residue = old+(1-old)*phase.real
    assert arb(".05") < transfer < arb(".053")
    assert abs(altered+(1-altered)*(phase+phase.conjugate())/2) < arb("1e-40")
    print("Selecting angle cells",count,"universal singlet real part > 0.005",
          "minimum required change of neutral population > 0.004")
    print("Seven triangles: real sum, protected neutral population and exact transfer",
          residue.str(16),altered.str(16),transfer.str(16))


if __name__=="__main__":
    verify()
