"""The heavy graph sees the arithmetic mean of both fitted positive endpoint energies."""
from pathlib import Path
from types import FunctionType
import runpy

from flint import acb, acb_mat, arb
import numpy as np

root = Path(__file__).resolve().parent
quarks = runpy.run_path(str(root/"ch39-one-tev-quark-fit.py"))
gauge = runpy.run_path(str(root/"ch39-common-height-gauge.py"))
metric,finite = quarks["metric"],quarks["finite"]
fit = metric["fit"]
shell = gauge["spectrum"]["shell"]
radius = arb("5e-8")
clocks = [v+arb(0,radius) for v in metric["common"]["angles"]]
model = fit["model"](*clocks)

def branch_weights(w):
    k = fit["evaluate"](model,w[:4],(0,0))[2]
    a = finite["fix_phase"](k,(7,6,8),arb.pi()-arb(".002"))
    b = finite["fix_phase"](k,(11,12,18),arb(".03606860"))
    return w+[a,b]

def source_bounds():
    hull = None
    # Chebyshev inputs x=cos(t) cover both closed interpolation pieces.
    # Bounding cos(j*t) directly keeps the interval evaluation well conditioned.
    for rows in (quarks["parts"][:22],quarks["parts"][22:]):
        for i in range(32):
            t = arb.pi()*(arb(i).union(arb(i+1)))/32
            terms = [(j*t).cos() for j in range(13)]
            w = [sum(a*b for a,b in zip(row,terms)) for row in rows]
            a,b = branch_weights(w[:11]),branch_weights(w[11:])
            mean = [(v+z)/2 for v,z in zip(a,b)]
            hull = mean if hull is None else [v.union(z) for v,z in zip(hull,mean)]
    assert min(hull)>0
    return hull

def edge_indices():
    pairs = {}
    for i,j,d,_ in shell["edges"]:
        if d in (724,804):
            pairs[i,j] = 2*int(i not in fit["even"])+int(d==804)
    for i,j in enumerate(fit["inner"]):
        for shift in range(3):
            k = fit["outer"][(i+shift)%3]
            pairs[tuple(sorted((j,k)))] = 10 if shift==0 else 4+3*(shift-1)+i
    for k,cycle in enumerate(((7,6,8),(11,12,18)),11):
        for i,j in zip(cycle,cycle[1:]+cycle[:1]):
            pairs[tuple(sorted((i,j)))] = k
    assert len(pairs)==27 and set(pairs)<=set((i,j) for i,j,_,_ in shell["edges"])
    return pairs

def static(ww):
    p,q = clocks
    c = gauge["spectrum"]["cos_axis"]
    v = gauge["spectrum"]["sin_axis"]
    su,cu = p.sin_cos()
    sv,cv = q.sin_cos()
    U = acb_mat([[acb(cu,su),0],[0,acb(cu,-su)]])
    V = acb_mat([[acb(cv,c*sv),v*sv],[-v*sv,acb(cv,-c*sv)]])
    I = acb_mat([[1,0],[0,1]])
    H = acb_mat(38,38)
    for i,w in enumerate(shell["V"]):
        for a in range(2):
            H[2*i+a,2*i+a] = w
    ids = edge_indices()
    for i,j,d,z in shell["edges"]:
        w = 1+ww[ids[i,j]] if (i,j) in ids else arb(1)
        G = U if d==724 else V if d==804 else I
        for a in range(2):
            for k in (i,j):
                H[2*k+a,2*k+a] += w
            for b in range(2):
                entry = -w*z*G[a,b]
                H[2*i+a,2*j+b] = entry
                H[2*j+b,2*i+a] = entry.conjugate()
    return H

def enclose_static(H):
    approximant = np.array([[complex(H[i,j].mid()) for j in range(38)] for i in range(38)])
    # Float eigenvalues propose cuts only. Inertia below every distinct cut is certified.
    centers = np.linalg.eigvalsh(approximant)[:18]
    intervals = []
    for i,v in enumerate(centers):
        lo,hi = arb(str(v))-arb(".001"),arb(str(v))+arb(".001")
        assert inertia(H,lo)<=i and inertia(H,hi)>i
        intervals.append(lo.union(hi))
    return intervals

def inertia(H,x):
    K = acb_mat(H)
    ids = list(range(K.nrows()))
    for i in ids:
        K[i,i]-=x
    total = 0
    while ids:
        j = max(ids,key=lambda k:float(abs(K[k,k].real).mid()))
        ids.remove(j)
        pivot = K[j,j].real
        assert pivot<0 or pivot>0
        total += int(pivot<0)
        for col,i in enumerate(ids):
            factor = K[i,j]/pivot
            for k in ids[col:]:
                v = K[i,k]-factor*K[j,k]
                K[i,k] = v.real if i==k else v
                if i!=k:
                    K[k,i] = v.conjugate()
    return total

def mass_bounds():
    weights = source_bounds()
    low = enclose_static(static([w.lower() for w in weights]))
    high = enclose_static(static([w.upper() for w in weights]))
    beta = gauge["beta"]
    limit,d = arb(200),gauge["data"]["d"]
    k = arb.pi()*(limit/d).sqrt()
    # tan(x)/x increases, so this bounds K(0)-sI-K(s) divided by s at every s<=200.
    slope = 18*k*(k/2).tan()/limit
    assert slope<arb(".00045") and max(high)<limit
    result = []
    for lo,hi in zip(low,high):
        cont = lo.lower()/(1+slope)
        result.append((cont/(1+beta*(cont+10))).lower().union(hi.upper()))
    return result

def setup():
    low = mass_bounds()
    r = gauge["r"]
    logs = [r.log()-x.log()/2 for x in low]
    assert sum(logs)<95 and max(logs)<6
    fn = gauge["integrate"]
    run = FunctionType(fn.__code__,dict(fn.__globals__,Ls=logs,WL=sum(logs)),argdefs=fn.__defaults__)
    moments = FunctionType(gauge["moments"].__code__,dict(gauge["moments"].__globals__,low=low))
    faces = FunctionType(gauge["faces"].__code__,dict(gauge["faces"].__globals__,exposures=logs,
                         left_sum=sum(logs),integrate=run,moments=moments))
    return low,faces

if __name__ == "__main__":
    low,faces = setup()
    print("extreme low modes",low[0],low[-1])
    faces("26.872","26.997","2.195","2.32")
