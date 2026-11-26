"""First variations at the 719th-zero window, with a closed-walk error bound."""
from pathlib import Path
import runpy
from flint import arb,arb_mat

root = Path(__file__).resolve().parent
sus = runpy.run_path(str(root/"ch39-neutral-clock-susceptibility.py"))
neutral = runpy.run_path(str(root/"ch39-su3-neutral-tower.py"))
q,r,raw = (sus[k] for k in ("q","r","raw"))
lift = neutral["lift"]
R = {d:arb_mat(a) for d,a in neutral["turn"].items()}
axes = (724,804)
def derivative(y,z,rot):
    cross = arb_mat([[0,-z,y],[z,0,0],[-y,0,0]])
    return -2*cross*rot
dR = [derivative(0,1,R[724]), derivative(lift["v"],lift["c"],R[804])]
zero = q(0)

def matrix_pair(zlo,zhi,upper):
    z = zlo if upper else zhi
    diag = [z for _ in range(19)]
    K = arb_mat(19,19)
    matrix = arb_mat(57,57) if not upper else None
    clock_edges = [[],[]]
    for i,j,d,ell,w in raw:
        if upper:
            y = (z.sqrt()*ell.lower()).lower()
            h = (y/y.tanh()+w.lower()).lower()
            b = (y/y.sinh()+w.upper()).upper()
        else:
            y = (z.sqrt()*ell.upper()).upper()
            h = (y/y.tanh()+w.upper()).upper()
            b = (y/y.sinh()+w.lower()).lower()
        diag[i]+=h
        diag[j]+=h
        K[i,j]=K[j,i]=-b
        rot = R.get(d,R[0])
        if matrix is not None:
            for a in range(3):
                for t in range(3):
                    matrix[3*i+a,3*j+t]=matrix[3*j+t,3*i+a]=-b*rot[a,t]
        if d in axes:
            clock_edges[axes.index(d)].append((i,j,b))
    for i,v in enumerate(diag):
        K[i,i]=v
        if matrix is not None:
            for a in range(3):
                matrix[3*i+a,3*i+a]=v
    G = K.inv()
    assert all(sum(G[i,j] for j in range(19))>0 for i in range(19))
    marks = [2*sum(b*(G[i,j]-b/(diag[i]*diag[j])) for i,j,b in edges) for edges in clock_edges]
    if upper:
        return marks
    twisted = matrix.inv()
    deriv = [-2*sum(b*dR[k][a,t]*twisted[3*i+a,3*j+t]
                    for i,j,b in edges for a in range(3) for t in range(3))
             for k,edges in enumerate(clock_edges)]
    return marks,deriv

def check():
    start = q("0.25")
    mult = 32*arb.pi()**2*r**4
    degree = sus["pot"]["short"]["base"]["arc_budget"]()[1].upper()/2
    weights = [sum(1+w.upper() for i,j,d,ell,w in raw if d==axis) for axis in axes]
    # Diamagnetic domination of every 3x3 Green block by its scalar counterpart.
    err0 = [8*w*(start/19+18*start**2/(19*16)) for w in weights]
    ratios = q(2)**(q(1)/64)
    lo = start
    out = [arb(0,v.upper()) for v in err0]
    for _ in range(64*30):
        hi = lo*ratios
        high = matrix_pair(lo,hi,True)
        low,vals = matrix_pair(lo,hi,False)
        width = (hi**2-lo**2)/2
        for i in range(2):
            # Each marked occurrence differentiates a rotation with trace norm four.
            err = 4*(high[i]-low[i]).upper()
            out[i] += width*(vals[i]+arb(0,err))
        lo = hi
    # 4 sum_n trace(B^n), n>=3 and ||B||<=degree/z.
    tail = 152*degree**3/lo
    vals = [(x+arb(0,tail.upper()))/mult for x in out]
    for k,v in zip(axes,vals):
        print(k,"clock energy derivative",v.lower().str(12),v.upper().str(12))
    theta = neutral["shell"]["theta"]
    a,b = (theta(d)-theta(244) for d in axes)
    slope = a*vals[0]+b*vals[1]
    print("arithmetic-height slope",slope.lower().str(12),slope.upper().str(12))
    assert vals[0]>0 and vals[1]<0 and slope<0

if __name__=="__main__":
    check()
