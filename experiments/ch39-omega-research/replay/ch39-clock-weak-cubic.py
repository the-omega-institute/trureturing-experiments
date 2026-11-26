"""Two neutral-clock Wilson factors in the weak loop masses: cubic spectrum and brightness."""
from itertools import permutations, product
from pathlib import Path
import runpy
from flint import acb,acb_mat,arb,fmpz_mat

root=Path(__file__).resolve().parent
cur=runpy.run_path(str(root/"ch39-weak-current-strength.py"))
pent=runpy.run_path(str(root/"ch39-five-clock-selection.py"))
q,short,raw=cur["q"],cur["short"],cur["raw"]
sites=cur["doublet"]["sites"]
rim=[i for i in range(19) if i!=14]
chords=[(i,j,d) for i,j,_,_,d in raw if 14 not in (i,j)]
P=fmpz_mat(len(chords),len(rim))
for n,(i,j,d) in enumerate(chords):
    P[n,rim.index(i)]=1
    P[n,rim.index(j)]=-1
outer=fmpz_mat(len(rim),1)
for perm in permutations(range(3)):
    v=sites.index(tuple(sites[0][k] for k in perm))
    outer[rim.index(v),0]=(-1)**sum(perm[i]>perm[j] for i in range(3) for j in range(i+1,3))
assert P.transpose()*P*outer==10*outer
for kind in (724,804):
    col=fmpz_mat([[outer[rim.index(i),0] if d==kind else 0] for i,j,d in chords])
    assert col.transpose()*col==fmpz_mat([[6]]) and P.transpose()*col==outer
assert all(d not in (724,804) for i,j,_,_,d in raw if 14 in (i,j))
s,t=[x.sin()**2 for x in pent["target"]]

def cubic(z):
    return ((z-13-2*(s+t))*z+23+13*(s+t)+3*s*t)*z-11*(1+s)*(1+t)

def roots():
    vals=[]
    for a,b in (("1.55","1.57"),("1.85","1.88"),("12.98","13.02")):
        lo,hi=q(a),q(b)
        sign=1 if cubic(hi)>0 else -1
        assert sign*cubic(lo)<0<sign*cubic(hi)
        for _ in range(36):
            mid=(lo+hi)/2
            if sign*cubic(mid)>0:
                hi=mid
            else:
                assert sign*cubic(mid)<0
                lo=mid
        vals.append(lo.union(hi))
    return vals

# Aggregate exact cyclic character factors by chord length. Direct matrix
# interval evaluation then preserves cancellations among symmetric edges.
omega=cur["omega"]
assign=cur["assignment"]
blocks={}
clock_lengths={}
for i,j,l,band,kind in raw:
    if kind not in blocks:
        blocks[kind]=(l,band,acb_mat(6,6),acb_mat(6,6))
    length,group,D,A=blocks[kind]
    for v in (i,j):
        if v in assign:
            a,_=assign[v]
            D[a,a]+=1
    if i in assign and j in assign:
        a,p=assign[i]
        b,v=assign[j]
        c=p.conjugate()*v
        A[a,b]+=c
        A[b,a]+=c.conjugate()

def ratio_bound(x,f):
    H=acb_mat(6,6)
    off={}
    for kind,(length,group,D,A) in blocks.items():
        y=x*length*f[group]/2
        sn,cs=y.sin_cos()
        off[kind]=y/sn
        H+=D*(off[kind]*cs)-A*off[kind]
    small=acb_mat([[H[a,b] for b in range(1,6)] for a in range(1,6)])
    rest=small.solve(acb_mat([[-H[a,0]] for a in range(1,6)]))
    inner=rest[1,0]
    assert cur["reps"][2]==6
    # Keep the common inner amplitude out of the quotient.
    phase=inner.imag/inner.real
    u,v=assign[7][1],assign[8][1]
    return off[804]/off[724]*(v.imag+v.real*phase)/(u.imag+u.real*phase)

def split(v,n):
    l,h=v.lower(),v.upper()
    return [l+(h-l)*q(k).union(q(k+1))/n for k in range(n)]

def check():
    rr=roots()
    dark=[-s*z/((z-1)*(z-11)-s*(2*z-11)) for z in rr]
    limits=[d.lower().union(d.upper()) for d in dark]
    results=[]
    base=[*short["fsbox"][:3],q("1.19576").union(q("1.21744"))]
    stack=[(base,0)]
    while stack:
        box,depth=stack.pop()
        f=short["lengths"](box[:3])
        if not f[3].overlaps(short["fsbox"][3]):
            continue
        f[3]=f[3].intersection(short["fsbox"][3])
        val=ratio_bound(box[3],f)
        if val>q("-1.19") and val<q("-.73"):
            results.append(val)
            continue
        assert depth<32, val
        # The three nontrivial widths are the two large optical bands plus mass.
        idx=(1,2,3,0)[depth%4]
        for seg in split(box[idx],2):
            child=box[:]
            child[idx]=seg
            stack.append((child,depth+1))
    ratio=min(z.lower() for z in results).union(max(z.upper() for z in results))
    assert all(not ratio.overlaps(d) for d in limits)
    pmax,qmax=q(".19955"),q(".21555")
    rest=1-9*(pmax+qmax)/8-(pmax*qmax).sqrt()/4
    assert rest>q(".48")
    print("mass squares / Q^2",*[z.str(12) for z in rr])
    print("dark ratios",*[d.str(10) for d in dark])
    print("current ratio",ratio.lower().str(12),ratio.upper().str(12),"boxes",len(results))

if __name__=="__main__":
    check()
