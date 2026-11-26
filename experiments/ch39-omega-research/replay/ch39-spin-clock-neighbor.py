"""Changes of charged fermion and boson determinants at neighboring arithmetic heights."""
from itertools import combinations, permutations
from pathlib import Path
import runpy
from flint import acb,acb_mat,arb,arb_mat

root = Path(__file__).resolve().parent
m = runpy.run_path(str(root/"ch39-su3-root-clock-lift.py"))
short,sh,q = m["short"],m["sh"],m["q"]
axes = (724,804)
ref = 719
def rotations(n):
    theta = sh["theta"]
    h = acb.zeta_zero(n).imag
    angles = [h*(theta(d)-theta(244)) for d in axes]
    (s,c),(t,b) = [v.sin_cos() for v in angles]
    v,a = m["v"],m["c"]
    return {724:acb_mat([[acb(c,s),0],[0,acb(c,-s)]]),
            804:acb_mat([[acb(b,a*t),v*t],[-v*t,acb(b,-a*t)]])}, angles

def rows(side):
    fs = short["fsbox"]
    V = sh["V"][:]
    if side=="right":
        for group,p in zip(short["m"]["groups"],short["shiftbox"]):
            for i in group:
                V[i]+=p
    raw = []
    for i,j,d,z in sh["edges"]:
        w = short["base"]["weights"][short["m"]["source_ids"][i,j]] if (i,j) in short["m"]["source_ids"] else q(0)
        ell = (sh["rho"]*sh["theta"](d)).sqrt()
        if side=="right":
            ell *= fs[short["band"](d)]
        raw.append((i,j,d,z,ell,w))
    return raw,V

def bounds(z,upper,raw,V):
    diag = [z+(v.lower() if upper else v.upper()) for v in V]
    K = arb_mat(19,19)
    bs = {}
    for i,j,d,phase,ell,w in raw:
        if z.is_zero():
            dy=by=q(1)
        else:
            y = (z.sqrt()*ell.lower()).lower() if upper else (z.sqrt()*ell.upper()).upper()
            dy,by=y/y.tanh(),y/y.sinh()
        if upper:
            h,b = (dy+w.lower()).lower(),(by+w.upper()).upper()
        else:
            h,b = (dy+w.upper()).upper(),(by+w.lower()).lower()
        diag[i]+=h
        diag[j]+=h
        bs[i,j]=b
        K[i,j]=K[j,i]=-b
    for i,v in enumerate(diag):
        K[i,i]=v
    G=K.inv()
    assert all(sum(G[i,j] for j in range(19))>0 for i in range(19))
    step = arb_mat(19,19)
    for i,j,d,phase,ell,w in raw:
        step[i,j] = bs[i,j]/diag[i]
        step[j,i] = bs[i,j]/diag[j]
    square = step*step
    small = step+square+square*step
    marks = [2*sum(bs[i,j]*(G[i,j]-small[i,j]/diag[j]) for i,j,d,_,_,_ in raw if d==k)
             for k in axes]
    return diag,bs,marks,step

def log_ratio(diag,bs,raw,U0,U):
    H = acb_mat(38,38)
    for i,x in enumerate(diag):
        H[2*i,2*i]=H[2*i+1,2*i+1]=x
    for i,j,d,phase,_,_ in raw:
        v = -phase*bs[i,j]
        if d not in axes:
            for k in (0,1):
                H[2*i+k,2*j+k]=v
                H[2*j+k,2*i+k]=v.conjugate()
    det=[]
    for rotations in (U0,U):
        for i,j,d,phase,_,_ in raw:
            if d in axes:
                mat=rotations[d]
                for a in range(2):
                    for b in range(2):
                        v=-phase*bs[i,j]*mat[a,b]
                        H[2*i+a,2*j+b]=v
                        H[2*j+b,2*i+a]=v.conjugate()
        h=H.det()
        assert h.imag.contains(0) and h.real>0
        det.append(h.real)
    return (det[1]/det[0]).log()

def short_cycles(raw,U0,U):
    links = {(i,j):(d,z) for i,j,d,z,_,_ in raw}
    coeffs = []
    ident = m["I"]
    for n in (3,4):
        for group in combinations(range(19),n):
            for perm in permutations(group[1:]):
                if perm[0]>perm[-1]:
                    continue
                cycle = [group[0],*perm]
                directed = tuple(zip(cycle,cycle[1:]+cycle[:1]))
                pairs = [tuple(sorted(edge)) for edge in directed]
                if not all(pair in links for pair in pairs) or not any(links[pair][0] in axes for pair in pairs):
                    continue
                products = []
                for turns in (U0,U):
                    rot = ident
                    for pair,edge in zip(pairs,directed):
                        d,z = links[pair]
                        u = turns.get(d,ident)*z
                        rot = rot*(u if edge==pair else u.conjugate().transpose())
                    products.append((rot[0,0]+rot[1,1]).real)
                val = 2*(products[0]-products[1])
                coeffs.append((directed,val))
    return coeffs

def interval(lo,hi,raw,V,U0,U,delta,cycles):
    gd,weight,low,sl=bounds(hi,False,raw,V)
    _,_,top,su=bounds(lo,True,raw,V)
    val = log_ratio(gd,weight,raw,U0,U)
    plus = minus = q(0)
    for cycle,c in cycles:
        left=right=q(1)
        for i,j in cycle:
            left*=sl[i,j]
            right*=su[i,j]
        val_inc = (right-left)*c
        plus += max(q(0),val_inc.upper())
        minus += min(q(0),val_inc.lower())
    err = sum(d*(h-l) for d,h,l in zip(delta,top,low)).upper()
    return (hi*hi-lo*lo)/2*(val+minus.lower().union(plus.upper())+arb(0,err))

def cert(side,n,steps):
    raw,V=rows(side)
    U0,a0=rotations(ref)
    U,a=rotations(n)
    ds = [4*abs(((x-y)/2).sin()) for x,y in zip(a,a0)]
    cycles = short_cycles(raw,U0,U)
    lo = q(".0625")
    out = interval(q(0),lo,raw,V,U0,U,ds,cycles)
    ratio=q(2)**(q(1)/steps)
    for _ in range(steps*35):
        hi=lo*ratio
        out+=interval(lo,hi,raw,V,U0,U,ds,cycles)
        lo=hi
    degree=short["base"]["arc_budget"]()[1].upper()/2
    # At large z all diagonal potentials >=-10, with 2*degree/z < 1/2.
    assert lo>4*degree+20
    tail=4*19*2*(2*degree)**3/lo
    return (out+arb(0,tail.upper()))/(32*arb.pi()**2*short["r"]**4)

def check():
    for n in (718,720):
        L=cert("left",n,16)
        R=cert("right",n,16)
        whole = 16*(R-L)
        assert whole < (q("2.49872e-12") if n==718 else q("1.32837e-12"))
        for tag,value in (("left",L),("right",R),("16(R-L)",whole)):
            print(n,tag,value.lower().str(12),value.upper().str(12),flush=True)

if __name__=="__main__":
    check()
