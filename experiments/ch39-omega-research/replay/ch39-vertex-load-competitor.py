"""A lower Gaussian action with the same vertex-load multisets and both direction budgets."""
from pathlib import Path
import runpy

from flint import arb

root = Path(__file__).resolve().parent
base = runpy.run_path(str(root/"ch39-one-tev-quark-fit.py"))
metric,finite = base["metric"],base["scale"]["finite"]
fit = metric["fit"]
pairs = [(i,fit["outer"][(j+k)%3]) for k in (1,2) for j,i in enumerate(fit["inner"])]

def matrix(w):
    k = fit["evaluate"](metric["model"],w[:4],(0,0))[2]
    for (i,j),val in zip(pairs,w[4:10]):
        finite["add_edge"](k,i,j,val)
    for i,j in zip(fit["inner"],fit["outer"]):
        finite["add_edge"](k,i,j,w[10])
    finite["fix_phase"](k,(7,6,8),arb.pi()-arb(".002"))
    finite["fix_phase"](k,(11,12,18),arb(".03606860"))
    return k

def exchange(w):
    a,b = w[4:7],w[7:10]
    inner = [a[i]+b[i] for i in (1,0,2)]
    outer = [a[(j-1)%3]+b[(j-2)%3] for j in range(3)]
    d1,d2 = outer[2]-inner[0], outer[0]-inner[1]
    a0 = (sum(a)-2*d1-d2)/3
    first = [a0,a0+d1,a0+d1+d2]
    return w[:4]+first+[v-z for v,z in zip(inner,first)]+w[10:]

def check():
    ratios = []
    minima = []
    for i in range(metric["ncell"]):
        q = metric["qlo"]+(arb(i)+arb(1)/2)*metric["step"]
        w = base["weights"](q)[:11]
        alt = exchange(w)
        assert min(alt)>0
        k = matrix(alt)
        deg,margin = finite["flux"]["degrees"](k)
        assert deg == (3,[3]) and margin>0
        ratio = matrix(w).det().real/k.det().real
        assert ratio>arb("1.0000034")
        ratios.append(ratio.lower())
        minima.append(min(alt[4:10]).lower())
    print("determinant ratio and alternate arc stiffness lower",min(ratios).str(12),min(minima).str(10))

if __name__ == "__main__":
    check()
