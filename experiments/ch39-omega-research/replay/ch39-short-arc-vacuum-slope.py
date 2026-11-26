"""Variation of the finite volume at fixed s and the two intermediate arc lengths."""
from pathlib import Path
import runpy
from flint import acb,arb_mat

m = runpy.run_path(str(Path(__file__).with_name("ch39-short-arc-vacuum.py")))
q,r,w = m["q"],m["r"],m["w"]
# All quantities in the derivative calculation use Q as the mass unit.
dprime = [row[0]/r for row in m["grad"]]
fprime = [q(1),q(0),q(0),-m["ratio"][0]]
a = [h/r**2 for h in m["mean_A"]]
n = m["count"]
delta,T = m["delta"],m["T"]

def mass_derivative(d,p):
    h = [v*v-b-z for v,b,z in zip(d,a,p)]
    A = arb_mat([[n[i]*row[i] for i in (2,3,4)] for row in ([q(1)]*5,d,h)])
    rhs = arb_mat([[sum(N*v*f for N,v,f in zip(n,row,dprime))] for row in
                   (d,h,[v*(v*v-2*(b+z)) for v,b,z in zip(d,a,p)])])
    if A.det().contains(0):
        return None
    # Eliminate by Cramer's rule. The derivative is unique on the positive-root branch.
    den = A.det()
    deriv = [q(0),q(0)]
    for col in range(3):
        B = arb_mat(A)
        for i in range(3):
            B[i,col] = rhs[i,0]
        deriv.append(B.det()/den)
    return deriv

def vertex_derivative(d,p,dp):
    # The rational term drops out after summing over both roots and all vertices,
    # by differentiating the four cancelled polynomial moments.
    total = q(0)
    for j,g in enumerate(m["m"]["groups"]):
        f = d[j]
        for i in g:
            b = m["vac"]["a"][i]+p[j]
            im = (b-f*f/4).sqrt()
            root = acb(f/2,im)
            dr = acb(dprime[j]/2,(dp[j]-f*dprime[j]/2)/(2*im))
            total += 8*(root**3*dr*(T+root).log()).real
    return total

def estimate(parts):
    x = m["coords"](parts[:3],parts[3])
    fs,inc = m["lengths"](x),m["increments"](x)
    p = [z/r**2 for z in m["branch"](inc,parts[3])]
    d = [(a+b)/r for a,b in zip(m["length"],inc)]
    dp = mass_derivative(d,p)
    if dp is None:
        return None
    V = vertex_derivative(d,p,dp)
    arcs = 6*m["zeta5"]*sum(k*v/f**5 for k,v,f in zip(m["inv4"],fprime,fs))
    pa = max(abs(v).upper() for v in dp)
    pd = max(abs(v).upper() for v in dprime)
    sf = max(abs(v/f).upper() for v,f in zip(fprime,fs))
    ellmax2 = q("1.018")**2*m["shell"]["rho"]*m["shell"]["theta"](804)
    low = 38*delta*(pa+18*sf*ellmax2*delta)
    off = m["hs"]*(pa/delta+2*pd/T)
    ip = 28*q(3).zeta()/r**2
    arc_change = 4*ip*sum(k*abs(v/f)/f**2 for k,v,f in zip(m["inv2"],fprime,fs))
    arc_inverse = ip*sum(k/f**2 for k,f in zip(m["inv2"],fs))*(pa/delta+pd/T)
    err = (low+off+arc_change+arc_inverse).upper()
    return (V+arcs-err).lower().union((V+arcs+err).upper())

def check():
    parts = [m["full"]]*3+[w["s0"].union(w["s1"])]
    todo = [(parts,0)]
    verified = []
    while todo:
        parts,depth = todo.pop()
        val = estimate(parts)
        if val is not None and val>5 and val<55:
            verified.append(val)
            continue
        assert depth<18, (depth,val)
        col = depth%4
        ax = parts[col]
        mid = ax.mid()
        for seg in (ax.lower().union(mid),mid.union(ax.upper())):
            sub = parts[:]
            sub[col] = seg
            todo.append((sub,depth+1))
    lo = min(v.lower() for v in verified)
    hi = max(v.upper() for v in verified)
    print("volume slope with respect to shortest length",lo.str(12),hi.str(12),"boxes",len(verified))

if __name__=="__main__":
    check()
