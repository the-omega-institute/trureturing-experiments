"""固定四个补空间见证时，六条低线的共同约束使第七壳真空读数偏离目标。"""
from pathlib import Path
from runpy import run_path
from flint import acb, acb_mat, arb

root = Path(__file__).resolve().parent
def lib(name):
    return run_path(root/("ch39-"+name+".py"))
witnesses = lib("six-common-spectrum-moment-witnesses")["WITNESSES"]
sieve = lib("six-cyclic-line-sieve")
adj, sq, row = [sieve[key] for key in ("adj", "squared", "row")]
heat = lib("chiral-minimal-neutral-gravity")["expint"]
def check(root, levels, blocks, model, frames):
    a0 = levels[0]
    trace = sq(blocks[0])
    T = sum(sq(m) for m in blocks[1:])
    fields = lib("three-index-free-atoms")["broken_fields"](root,range(2306888100092,2306888100098))[2]
    eps = (arb(3000)/arb("1.22089e19")-a0).sqrt()
    print("trace, T, det", trace.str(10),T.str(10),(adj(blocks[0])*blocks[0]).det().str(10),flush=True)
    def gravity(lo,hi,field):
        k = w = arb(0)
        for g,a,z in field:
            x,y = lo**2,hi**2
            k += a*(x*heat(2,z/x)).union(y*heat(2,z/y))/(2*arb.pi())
            w -= g*(x*x*heat(3,z/x)).union(y*y*heat(3,z/y))/(4*arb.pi())
        return k,w
    windows = (("-.050","-.030"),(".034",".038"),(".030",".050"),("-.0121","-.0117"))
    for (labels,re,im),(left,right) in zip(witnesses,windows,strict=True):
        B = acb_mat([row(s["zero"],i-1).tolist()[0] for s,i in zip(model,labels)])
        J = adj(B)*(B*adj(B)).inv()
        F = acb_mat([[int(i==j) for j in range(7)] for i in range(7)]) - J*B
        k = max(range(7),key=lambda j:F[j,j].real.lower())
        n = acb_mat([[F[j,k]] for j in range(7)])
        n = n/sq(n).sqrt()
        cap = eps*sum(sq(acb_mat([[J[i,j]] for i in range(7)])).sqrt() for j in range(6))
        print("LABELS, CAP",labels,cap.str(10),flush=True)
        q=acb_mat([[acb(x,y)] for x,y in zip(re,im)])
        norm2 = sum(x*x+y*y for x,y in zip(re,im))
        assert trace > cap**2
        for j,(s,f,field) in enumerate(zip(model[:1],frames[:1],fields[:1],strict=True)):
            K=W=arb(0)
            for i,w in enumerate(s["w"]):
                if i in s["triv"]:
                    r=row(s["zero"],i)
                    h=abs((r*n)[0,0])
                    err=eps*sum(abs(x) for x in (r*J).entries())
                    lo=a0+max(arb(0),(h*(trace-cap**2).sqrt()-err).lower())**2
                    hi=a0+(h*trace.sqrt()+err)**2
                else:
                    lo=hi=a0+T*abs((acb_mat([f[i]])*q)[0,0])**2/norm2
                u,v=gravity(lo.lower(),hi.upper(),field)
                K+=w*u
                W+=w*v
            assert arb(left) < W < arb(right)
            print(j+7,"K-1", (K-1).str(10),"W",W.str(10),flush=True)
source = lib("fifth-full-positive-joint")
lib("six-trivial-shared-endpoints")["experiment"](root,*source["verify"](root,return_operator="spectral"),check=check)
