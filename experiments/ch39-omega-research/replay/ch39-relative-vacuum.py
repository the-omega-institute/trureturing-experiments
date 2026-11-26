"""A continuous branch cancelling three endpoint coefficients and two finite-scale moments."""
from pathlib import Path
import runpy

from flint import arb, arb_mat

root = Path(__file__).resolve().parent
w = runpy.run_path(str(root / "ch39-weyl-balanced-arcs.py"))
base, shell, q, r = w["b"], w["shell"], w["q"], w["r"]
V = shell["V"]
types = [tuple(sorted(p)) for p in shell["sites"]]
orbits = list(dict.fromkeys(types))
groups = [[i for i,t in enumerate(types) if t == group] for group in orbits]
count = [len(group) for group in groups]
assert orbits == [(1,5,29), (13,13,23), (7,17,23), (11,11,25), (17,17,17)]
assert count == [6,3,6,3,1]

# The length and both of its derivatives are constant on each permutation orbit.
edge_classes = [[] for _ in V]
added = [q(0) for _ in V]
source_ids = w["source"]["edge_indices"]()
for i,j,d,_ in shell["edges"]:
    for k in (i,j):
        edge_classes[k].append(d)
        if (i,j) in source_ids:
            added[k] += base["weights"][source_ids[i,j]]
length, du, dv, mean_A, potential = [], [], [], [], []
for group in groups:
    ds = sorted(edge_classes[group[0]])
    # V is a function of the coordinate sum, hence constant on each orbit.
    assert all(sorted(edge_classes[i]) == ds for i in group)
    bands = [sum((shell["rho"]*shell["theta"](d)).sqrt() for d in ds if w["band"](d)==j) for j in range(3)]
    length.append(sum(bands))
    du.append(bands[0]-w["a0"]*bands[2])
    dv.append(bands[1]-w["a1"]*bands[2])
    potential.append(V[group[0]])
    mean_A.append(sum(V[i]+added[i] for i in group)/len(group))

def branch(u, v, s):
    delta = [a*(u-1)+b*(v-1) for a,b in zip(du,dv)]
    D = [d+e for d,e in zip(length,delta)]
    delta2 = [e*(f+d) for e,f,d in zip(delta,D,length)]
    mu = r**2*(-2*s).exp()
    mass_sum = sum(n*g/2 for n,g in zip(count,delta2)) - 9*mu
    # Subtract the common reference in tr((D_R-D_L)A) using exact length balance.
    ref = potential[2]
    d_weight = sum(n*e*((f*f+f*d+d*d)/3-(a-ref)) for n,e,f,d,a in zip(count,delta,D,length,mean_A))
    d_weight -= mu*(6*D[0]+3*D[1])
    gap = D[2]-D[3]
    assert gap > 0
    h0 = (d_weight-D[3]*mass_sum)/gap
    h1 = (D[3]-D[4])/gap
    # p = intercept + slope*z. The linear coefficients above are for orbit totals.
    intercept = [mu,mu,h0/6,(mass_sum-h0)/3,q(0)]
    slope = [q(0),q(0),h1/6,(-1-h1)/3,q(1)]
    qa = qb = qc = q(0)
    for n,d,f,e2,a,x,b in zip(count,length,D,delta2,mean_A,intercept,slope):
        linear = f*f-(a-ref)
        qa -= n*b*b/2
        qb += n*b*(linear-x)
        qc += n*(e2*(a-ref/2-(f*f+d*d)/4)+linear*x-x*x/2)
    disc = qb*qb-4*qa*qc
    assert disc > 0 and qb < 0 and qc > 0 and qa < 0
    z = 2*qc/(-qb+disc.sqrt())
    return [a+b*z for a,b in zip(intercept,slope)]

def endpoints(shifts):
    value = first = q(0)
    for n,V0,p in zip(count,potential,shifts):
        x, y = V0/r**2,p/r**2
        mult = 2*n*(-x).exp()
        change = (-y).expm1()
        value += mult*change
        first += mult*(x*change+y*(-y).exp())
    return value, first

def gap_check(shifts):
    eta = q("1.018")**2
    beta = base["data"]["beta"]*eta
    energy = r*r*(1+10*beta)/(1-r*r*beta)
    k = arb.pi()*(energy*eta/base["data"]["d"]).sqrt()
    high = min((v+p).lower() for v,p in zip(potential[2:],shifts[2:]))
    assert k<arb.pi() and high-energy-36*k/k.sin()>0 and high>shifts[0].upper()
    return high

def fit():
    center, radii = w["center"], w["radii"]
    s0,s1 = w["s0"], w["s1"]
    shifts = branch(*w["box"],s0.union(s1))
    print("shift bounds",[(float(p.lower()),float(p.upper())) for p in shifts[2:]])
    print("high onsite floor",gap_check(shifts).str(12))
    trace,_ = base["arc_budget"]()
    # The static graph stiffness is positive for both bundles. On x >= -0.00007,
    # -1.001 < (exp(-x))' <= 0 and -0.136 < (x exp(-x))' < 1.001.
    assert min(v.lower() for v in V)/r**2 > q("-.00007")
    errors = [(1+q("1.018")**2)*base[name]+w["tail"]+q(factor)*trace
              for name,factor in (("errS","1.001"),("errM","1.137"))]
    J = w["jac"](*center)
    P = arb_mat([[J[i,j].mid() for j in range(2)] for i in range(2)]).inv()
    for i in range(2):
        j = 1-i
        for sign in (-1,1):
            margins = []
            for k in range(32):
                f = center[:]
                f[i] += sign*radii[i]
                f[j] += radii[j]*(q(k).union(q(k+1))/16-1)
                interior = w["dirichlet"](*f)
                for part in range(12):
                    s = s0+(s1-s0)*q(part).union(q(part+1))/12
                    vertex = endpoints(branch(*f,s))
                    result = arb_mat([[d+h+arb(0,e.upper())] for d,h,e in zip(interior,vertex,errors)])
                    margins.append((sign*(P*result)[i,0]).lower())
            print("oriented thermal face",i,sign,min(margins).str(12))
            assert min(margins)>0
    w["gauge"]()

if __name__ == "__main__":
    fit()
