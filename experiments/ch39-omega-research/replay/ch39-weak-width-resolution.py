"""A uniform partial half-width in squared-mass units for the highest clock-sensitive line."""
from pathlib import Path
import runpy
from flint import acb,acb_mat,arb

root=Path(__file__).resolve().parent
labels=runpy.run_path(str(root/"ch39-weak-clock-labels.py"))
width=runpy.run_path(str(root/"ch39-weak-clock-width.py"))
m,q=labels["m"],labels["q"]
short=m["short"]
z=width["z3"]
r=m["s"]*z/((z-1)*(z-11)-m["s"]*(2*z-11))
norm=q(".9")*(r*r+1)-q(".2")*r+q("1.1")*(r+1)**2*((z-1)/(z-11))**2
mu=width["wv"]
phase=(1-4*mu**2/z).sqrt()*(2+((z-2*mu**2)/(2*mu**2))**2)/(8*arb.pi())
inv=width["upper_inverse"](z.sqrt())
gain=4*arb.pi()*phase/(inv*norm)
lo_mass=q("1.19576")
hi_mass=q("1.21744")
base=[*short["fsbox"][:3],lo_mass.union(hi_mass)]
def bound(box):
    f=short["lengths"](box[:3])
    if not f[3].overlaps(short["fsbox"][3]):
        return None
    f[3]=f[3].intersection(short["fsbox"][3])
    H=acb_mat(6,6)
    arcs={}
    offs={}
    for kind,(ell,g,D,A) in m["blocks"].items():
        L=ell*f[g]
        y=box[3]*L/2
        sn,cs=y.sin_cos()
        off=y/sn
        H+=D*(off*cs)-A*off
        arcs[kind]=((L/2)**2,(1-sn*cs/y)/(2*sn**2),(sn/y-cs)/(2*sn**2))
        offs[kind]=off
    minor=acb_mat([[H[i,j] for j in range(1,6)] for i in range(1,6)])
    rest=minor.solve(acb_mat([[-H[i,0]] for i in range(1,6)]))
    amp=[acb(1)]+[rest[i,0] for i in range(5)]
    u,v=labels["u"],labels["v"]
    inner=amp[2]
    # Form the linear combination in a shared endpoint amplitude; do not
    # subtract intervals for the final spectral weights of the two edges.
    diff=((r*offs[724]*u.imag+offs[804]*v.imag)*inner.real+
          (r*offs[724]*u.real+offs[804]*v.real)*inner.imag)
    jA=offs[724]*(u.imag*inner.real+u.real*inner.imag)
    by_ratio=jA.abs_lower()*(r.lower()-q("1.09974"))
    floor=max(diff.abs_lower(),by_ratio.lower())
    signed=m["assign"]
    vals=[acb(0)]*19
    for i,(a,p) in signed.items():
        vals[i]=amp[a]*p
    G=W=q(0)
    for i,j,_,_,kind in m["raw"]:
        w,d,c=arcs[kind]
        G+=(w*(d*(abs(vals[i])**2+abs(vals[j])**2)+2*c*(vals[i].conjugate()*vals[j]).real)).upper()
        W+=w.lower()
    return (gain.lower()*24*W*floor**2/G**2).lower()

def check():
    floor=q("1.0001e-7")
    stack=[(base,0)]
    vals=[]
    while stack:
        box,depth=stack.pop()
        v=bound(box)
        if v is None:
            continue
        if v>floor:
            vals.append(v)
            continue
        assert depth<30, (depth,v)
        idx=(1,2,3,0)[depth%4]
        for part in m["split"](box[idx],2):
            child=box[:]
            child[idx]=part
            stack.append((child,depth+1))
    print("uniform lower M Gamma_fg / Q^2",min(vals).str(12),"cells",len(vals))
if __name__=="__main__":
    check()
