"""在任意四站复三平面上，以投影交错检验共同的单位区间二次有理振幅。"""
from pathlib import Path
from runpy import run_path

from flint import acb, arb, ctx


MASS = (((".00000601",".00000629"),(".00306",".00316"),(".8573",".8659")),
        ((".0000133",".0000137"),(".000265",".000271"),(".01392",".01410")))


def spectrum(data,weights):
    p,t,re,im,plus,minus = [sum(v*w for v,w in zip(row,weights,strict=True)) for row in data]
    mid = (p+t)/2
    radius = abs(acb((p-t)/2,abs(acb(re,im))))
    # 分别对四个真特征值的上下端排序，以免在区间交叠时选择错误次序。
    eig = (mid-radius,mid+radius,plus,minus)
    return sorted(x.lower() for x in eig),sorted(x.upper() for x in eig)


def cell(theta,width,anchors,angles,distances,blocks,masses):
    norm_bounds = []
    for angle in angles:
        a = abs((angle-theta).sin())
        # ||v||=sqrt(sin²+η²cos²) 随 |sin|、η 分别非降。
        lo = max(arb(0),a.lower())
        hi = min(arb(1),a.upper())
        e0 = max(arb(0),width.lower())
        e1 = min(arb(1),width.upper())
        low = (lo*lo+(1-lo*lo)*e0*e0).sqrt().lower()
        high = (hi*hi+(1-hi*hi)*e1*e1).sqrt().upper()
        norm_bounds.append((low,high))
    keep = []
    for j in anchors:
        eps = (masses[j//19][0][1]/(masses[0][2][0]-masses[0][0][1])).sqrt().upper()
        cos = (1-eps*eps).sqrt().lower()
        lower,upper = [],[]
        for i,(d0,d1) in enumerate(norm_bounds):
            if i==j:
                lower.append(arb(0))
                upper.append(masses[j//19][0][1])
                continue
            j0,j1 = norm_bounds[j]
            den = d0*j0
            s0 = (max(arb(0),width.lower())*distances[i][j]/(d1*j1)).lower()
            s1 = min(arb(1),(width.upper()*distances[i][j]/den).upper()) if den>0 else arb(1)
            f0 = max(arb(0),(s0*cos-eps).lower())
            lower.append(((masses[0][2][0]-masses[0][0][1])*f0*f0).lower())
            upper.append(min(arb(1),(masses[0][0][1]+(s1+eps)**2).upper()))
        survives = True
        for q,B in enumerate(blocks):
            low = spectrum(B,lower[19*q:19*q+19])[0]
            if any(low[k]>masses[q][k][1] for k in range(3)):
                survives = False
                break
            high = spectrum(B,upper[19*q:19*q+19])[1]
            if any(high[k+1]<masses[q][k][0] for k in range(3)):
                survives = False
                break
        if survives:
            keep.append(j)
    return keep


def verify(root):
    native = run_path(str(root/"ch39-free-four-site-native-positive-spectrum.py"))
    shell = run_path(str(root/"ch39-three-ray-chiral-boundary.py"))
    ordered = run_path(str(root/"ch39-bounded-common-rational-amplitude.py"))["ordered"]
    ctx.prec = 192
    masses = [[[arb(v) for v in pair] for pair in row] for row in MASS]
    for n in (2_306_888_100_085,2_306_888_100_086):
        rows = native["native"](root,n)
        ctx.prec = 192
        h = acb.zeta_zero(n).imag
        angles = [(t-18).atan() for q in (-4,2) for t in ordered(shell,q,h)]
        distances = [[abs((a-b).sin()) for b in angles] for a in angles]
        blocks = []
        for qs in rows:
            coords = [(p[0,0].real,p[1,1].real,p[0,1].real,p[0,1].imag,
                       p[2,2].real+p[2,3].imag,p[2,2].real-p[2,3].imag) for p in qs]
            blocks.append(list(zip(*coords,strict=True)))
        work = [(0,0,0,list(range(38)))]
        cells,depth = 0,0
        while work:
            level,u,v,anchors = work.pop()
            step = arb(1)/2**level
            theta = ((u*step-arb(".5")).union((u+1)*step-arb(".5")))*arb.pi()
            width = (v*step).union((v+1)*step)
            keep = cell(theta,width,anchors,angles,distances,blocks,masses)
            cells += 1
            depth = max(depth,level)
            if any(j<19 for j in keep) and any(j>=19 for j in keep):
                assert level < 35,(n,level,u,v,keep)
                work.extend((level+1,2*u+i,2*v+j,keep) for i in range(2) for j in range(2))
        print("编号、分母形状盒数、最大二分深度",n,cells,depth,flush=True)


if __name__=="__main__":
    verify(Path(__file__).resolve().parent)
