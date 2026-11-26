"""区间覆盖原壳第三、第四本征值在五荷全钟面上的最小差。"""
from functools import cache
from pathlib import Path
from runpy import run_path
import sys

from flint import acb, acb_mat, arb, ctx


@cache
def phase(index, n):
    return acb(0, arb.pi()*(2*index+1)/n).exp()


def eigen_data(blocks, u, v, width, change):
    base, a, b = blocks
    a, b = a*u, b*v
    adj = lambda z:z.conjugate().transpose()
    h = base + a+adj(a)+b+adj(b)
    vals, vectors = h.eig(right=True, algorithm="rump")
    order = sorted(range(len(vals)), key=lambda j:float(vals[j].real))
    assert all(x.imag.contains(0) for x in vals)
    assert all(vals[i].real < vals[j].real for i,j in zip(order[:2],order[1:3]))
    out = []
    for col in order[:2]:
        lam = vals[col].real
        gap = min(abs(lam-vals[j].real).lower() for j in order if j!=col)
        vec = acb_mat([[vectors[i,col]] for i in range(vectors.nrows())])
        norm = (adj(vec)*vec)[0,0].real
        aa, bb = ((adj(vec)*c*vec)[0,0]/norm for c in (a,b))
        error = (change**2/(gap-2*change) if gap > 2*change else None)
        out.append((lam, aa, bb, error))
    return out


def verify(root):
    shell = run_path(str(root/"ch39-three-ray-chiral-boundary.py"))
    light = run_path(str(root/"ch39-native-light-projector-pump.py"))
    ctx.prec = 112
    adj = light["adjoint"]
    bases = list(light["rotation_bases"](shell))
    target = arb("1.61363")
    soft = 6*arb("1e-8")
    for q in (1,-4,2,-3,6):
        base, links = light["clock_matrices"](shell,q)
        blocks = [tuple(adj(b)*h*b for h in (base,*links)) for b in bases]
        if q == 6:
            star = [eigen_data(x, acb(0,arb("1.19205")).exp(), acb(0,arb("4.862")).exp(),arb(0),arb(0))
                    for x in blocks]
            high = min(row[1][0].upper() for row in star)
            low = max(row[0][0].lower() for row in star)
            upper = high-low
            assert upper<arb("1.613637")
        stack = [(i,j,8) for i in range(8) for j in range(8)]
        count,deepest = 0,8
        while stack:
            i,j,n = stack.pop()
            count += 1
            deepest = max(deepest,n)
            assert n <= 2**18
            width = arb.pi()/n
            change = 4*arb(2).sqrt()*(width/2).sin()
            at = [eigen_data(x,phase(i,n),phase(j,n),width,change) for x in blocks]
            safe = True
            for l in (x[0] for x in at):
                for h in (x[1] for x in at):
                    gap = h[0]-l[0]
                    if gap-2*change-soft > target:
                        continue
                    if h[3] is not None and l[3] is not None:
                        a,b = h[1]-l[1],h[2]-l[2]
                        drift = 2*width*(abs(a.imag)+abs(b.imag))+width**2*(abs(a)+abs(b))
                        if gap-drift-h[3]-l[3]-soft > target:
                            continue
                    safe = False
                    break
                if not safe:
                    break
            if not safe:
                stack.extend((2*i+di,2*j+dj,2*n) for di in (0,1) for dj in (0,1))
        print("荷、核验方格数、最细每轴分格数、全域谱隙严格下界",q,count,deepest,str(target),flush=True)
    print("荷六给出的全域最小谱隙上界",upper.str(13),flush=True)


if __name__=="__main__":
    root = Path(sys.argv[1]).resolve() if len(sys.argv)>1 else Path(__file__).resolve().parent
    del sys.argv[1:]
    verify(root)
