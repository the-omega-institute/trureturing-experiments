"""Cancel the finite relative one-loop volume by adjusting the three shortest arcs together."""
from functools import lru_cache
from itertools import product
from math import floor, ceil
from pathlib import Path
import runpy
from flint import arb, arb_mat

vac = runpy.run_path(str(Path(__file__).with_name("ch39-relative-vacuum-value.py")))
m,w,q,r,base = (vac["m"],vac["w"],vac["q"],vac["r"],vac["m"]["base"])
shell, count, potential, mean_A, length = (m[k] for k in ("shell","count","potential","mean_A","length"))
T,delta,a0,linear,slope = (vac[k] for k in ("T","delta","a0","linear","slope"))
prime,remainder,left_rem,left_d,a = (vac[k] for k in ("prime","remainder","left_rem","left_d","a"))

def band(d):
    return 0 if d==12 else 1+w["band"](d)

bands = [[q(0) for _ in range(4)] for _ in count]
total = [q(0) for _ in range(4)]
inv2,inv4 = total[:],total[:]
index = {i:j for j,g in enumerate(m["groups"]) for i in g}
for i,j,d,_ in shell["edges"]:
    ell = (shell["rho"]*shell["theta"](d)).sqrt()
    k = band(d)
    total[k] += ell
    inv2[k] += 1/(ell*r)**2
    inv4[k] += 1/(ell*r)**4
    for group in (index[i],index[j]):
        bands[group][k] += ell/count[group]
ratio = [b/total[3] for b in total[:3]]
grad = [[row[k]-ratio[k]*row[3] for k in range(3)] for row in bands]
counts = [sum(band(d)==k for _,_,d,_ in shell["edges"]) for k in range(4)]
assert counts[0] == 3

def lengths(x):
    return [*x,1-sum(c*(b-1) for c,b in zip(ratio,x))]

def increments(x):
    return [sum(c*(b-1) for c,b in zip(row,x)) for row in grad]

def branch(delta,s):
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

def vertex_integral(inc,shifts):
    p = [x/r**2 for x in shifts]
    result = q(0)
    for j,(g,d,pj) in enumerate(zip(m["groups"],left_d,p)):
        e = inc[j]/r
        f = d+e
        diff2 = e*(f+d)
        # All vertices of an orbit share d,p, but can have distinct fitted stiffnesses.
        for i in g:
            x = a[i]
            derivative = q(0)
            for k in range(4):
                b = x+pj*q(k).union(q(k+1))/4
                if j < 2:
                    derivative += (prime(b)-linear-slope*b)/4
                else:
                    axis = b.union(a0)+delta
                    diff = b-a0
                    derivative += diff*diff*(1/axis+delta/(axis*axis))/8
            # Restore the two subtracted polynomial moments using C2=C4=0.
            result += 2*(pj*derivative + linear*diff2/2
                        +slope*(diff2*(x-(f*f+d*d)/4)+f*f*pj)
                        +remainder(x+pj,f))
        result -= 2*left_rem[j]
    return result

center = list(map(q,(".948","1.006",".986")))
motion = list(map(q,(".087","-.010",".026")))
radii = list(map(q,(".020",".010",".008")))
s0,s1 = w["s0"],w["s1"]
def coords(t,s):
    return [c+d*(s-q("2.255"))+h*z for c,d,h,z in zip(center,motion,radii,t)]
full = q(-1).union(q(1))
domain = coords([full]*3,s0.union(s1))
fsbox = lengths(domain)
assert min(z.lower() for z in fsbox)>q(".91") and max(z.upper() for z in fsbox)<q("1.018")
incbox = increments(domain)
shiftbox = branch(incbox,s0.union(s1))
assert m["gap_check"](shiftbox)>q("800000")

# Uniform perturbation bound for the determinant integral, using the new shortest length.
c = 1-10/(delta*r*r)
arc = 28*q(3).zeta()*sum(k*(1+1/f.lower()**2) for k,f in zip(inv2,fsbox))/(r*r*c)
hs = 4*sum(base["weights"][idx]**2 for idx in m["source_ids"].values())/r**4
pa = max(abs(x).upper() for x in shiftbox)/r**2
pd = max(abs(x).upper() for x in incbox)/r
off = hs*(pa/delta+2*pd/T)/c**3
xmax = max((v+p).upper() for v,p in zip(potential,shiftbox))/r**2
_,norm = base["arc_budget"]()
assert xmax+delta+(norm+18*(q("1.018")*r*T*(shell["rho"]*shell["theta"](804)).sqrt()+2))/r**2<10
ev = (arc+off+38*delta**2*((10/delta).log()+q(".5"))/2).upper()
trace,_ = base["arc_budget"]()
es,em = [(1+q("1.018")**2)*base[name]+w["tail"]+q(err)*trace
         for name,err in (("errS","1.001"),("errM","1.137"))]
# Seven continuum harmonics suffice with the same previously bounded tail.
rows = [[] for _ in fsbox]
for d in w["classes"]:
    multiplicity = 2*sum(t==d for _,_,t,_ in shell["edges"])
    for j in range(1,8):
        k = (arb.pi()*j)**2/(r*r*shell["rho"]*shell["theta"](d))
        rows[band(d)].append((multiplicity,k))

P = [list(map(q,row)) for row in ((".105",".041",".041"),("-.094","-.064","-.00745"),
                                  (".00074",".052",".00216"))]
assert arb_mat(P).det()<0
zeta5 = q(5).zeta()

@lru_cache(maxsize=4096)
def arc_data(j,k0,k1):
    # The calculation uses exact integer grid boundaries.
    f = q(k0).union(q(k1))/10000000
    S=M=ds=dm=q(0)
    for n,k in rows[j]:
        inv = 1/(f*f)
        y = k*inv
        ch = (-(inv-1)*k).expm1()
        S += n*(-k).exp()*ch
        M += n*k*(-k).exp()*(inv*ch+inv-1)
        der = 2*n*y*(-y).exp()/f
        ds += der
        dm += der*(y-1)
    V = -q("1.5")*zeta5*inv4[j]*(1/f**4-1)
    dv = 6*zeta5*inv4[j]/f**5
    return (S,M,V),(ds,dm,dv)

def arc_residual(i,x):
    xm = [v.mid() for v in x]
    result = q(0)
    grad = []
    for j,(f,fm) in enumerate(zip(lengths(x),lengths(xm))):
        def enclosure(v):
            lo,hi = floor(float(v.lower())*10000000)-1,ceil(float(v.upper())*10000000)+1
            assert q(lo)/10000000<v<q(hi)/10000000
            return arc_data(j,lo,hi)
        vals = enclosure(fm)[0]
        ders = enclosure(f.union(fm))[1]
        result += sum(z*h for z,h in zip(vals,P[i]))
        grad.append(sum(z*h for z,h in zip(ders,P[i])))
    return result+sum((h-rat*grad[3])*(z-zm) for h,rat,z,zm in zip(grad[:3],ratio,x,xm))

def margin(i,sign,box,s):
    x = coords(box,s)
    inc = increments(x)
    shifts = branch(inc,s)
    S,M = m["endpoints"](shifts)
    V = vertex_integral(inc,shifts)
    err = sum(abs(v)*e for v,e in zip(P[i],(es,em,ev))).upper()
    value = arc_residual(i,x)+P[i][0]*S+P[i][1]*M+P[i][2]*V
    return sign*value+arb(0,err)

def check_face(i,sign):
    # Bisect only boxes whose direct oriented residual still contains zero.
    free = [j for j in range(3) if j!=i]
    pending = [([full if j!=i else q(sign) for j in range(3)],s0.union(s1),0)]
    margins = []
    while pending:
        box,s,depth = pending.pop()
        value = margin(i,sign,box,s)
        if value>abs(P[i][2])*q(".024"):
            margins.append(value.lower())
            continue
        assert depth<17, (i,sign,depth,value)
        k = depth%3
        axis = s if k==2 else box[free[k]]
        mid = axis.mid()
        for seg in (axis.lower().union(mid),mid.union(axis.upper())):
            child = box[:]
            if k<2:
                child[free[k]] = seg
            pending.append((child,seg if k==2 else s,depth+1))
    print("spectral face",i,sign,len(margins),min(margins).str(8),flush=True)

def check():
    for i,sign in product(range(3),(-1,1)):
        check_face(i,sign)
    # The gauge integrator uses the same uniform upper optical factor, low-site
    # floor and positive high block already checked for this new rectangle.
    w["gauge"]()
    print("length factors",[(float(f.lower()),float(f.upper())) for f in fsbox])

if __name__=="__main__":
    check()
