"""Enclose the one-loop Hessian and the small-energy zero density by positive walk comparison."""
from pathlib import Path
import runpy
from flint import arb,arb_mat,acb_mat

m = runpy.run_path(str(Path(__file__).with_name("ch39-neutral-clock-density.py")))
q,r,sh,pot = (m[k] for k in ("q","r","sh","pot"))
axes = (724,804)
raw = []
for i,j,d,_ in sh["edges"]:
    f = pot["short"]["fsbox"][pot["short"]["band"](d)]
    ell = (sh["rho"]*sh["theta"](d)).sqrt()*f
    raw.append((i,j,d,ell,pot["edges"][i,j][2]))

def response(zlo,zhi,upper):
    # All comparison quantities are point bounds, so inverting K creates negligible wrapping.
    diag = [zlo if upper else zhi for _ in range(19)]
    K = arb_mat(19,19)
    B = [arb_mat(19,19) for _ in axes]
    chords = [[],[]]
    for i,j,d,ell,w in raw:
        if upper:
            y = (zlo.sqrt()*ell.lower()).lower()
            diag_term = (y/y.tanh()+w.lower()).lower()
            b = (y/y.sinh()+w.upper()).upper()
        else:
            y = (zhi.sqrt()*ell.upper()).upper()
            diag_term = (y/y.tanh()+w.upper()).upper()
            b = (y/y.sinh()+w.lower()).lower()
        K[i,j]=K[j,i]=-b
        diag[i]+=diag_term
        diag[j]+=diag_term
        if d in axes:
            h = axes.index(d)
            B[h][i,j],B[h][j,i] = b,-b
            chords[h].append((i,j,b))
    for i,v in enumerate(diag):
        K[i,i]=v
    # K positive definite implies rho(D^-1/2 W D^-1/2)<1 by nonnegative symmetry.
    G = K.inv()
    assert all(sum(G[i,j] for j in range(19))>0 for i in range(19))
    T = [G*b for b in B]
    diag_h = [4*(2*sum(b*G[i,j] for i,j,b in c)+(t*t).trace()) for c,t in zip(chords,T)]
    mixed = 4*(T[0]*T[1]).trace()
    return [*diag_h,sum(diag_h)+2*mixed,sum(diag_h)-2*mixed]

def check():
    # The low-energy coercivity estimate follows from the positive scalar spectral gap.
    start = q(1)
    lap = arb_mat(19,19)
    for i,j,d,_,w in raw:
        v = 1+w.lower()
        lap[i,i]+=v
        lap[j,j]+=v
        lap[i,j]=lap[j,i]=-v
    # L + projector onto constants times 9 - 8 I is positive.
    test = arb_mat([[lap[i,j]+q(9)/19-8*int(i==j) for j in range(19)] for i in range(19)])
    assert pot["m"]["source"]["inertia"](acb_mat(test),q(0))==0
    weights = [sum((1+w.upper()) for i,j,d,_,w in raw if d==axis) for axis in axes]
    low_err = [8*w*(start/19+q(18)*start**2/(19*16)) for w in weights]
    lower = [q(0) for _ in range(4)]
    upper = [*low_err,sum(low_err),sum(low_err)]
    ratio = q(2).sqrt().sqrt().sqrt().sqrt()
    lo = start
    for _ in range(16*32):
        hi = ratio*lo
        width = (hi**2-lo**2)/2
        lower = [v+width*k for v,k in zip(lower,response(lo,hi,False))]
        upper = [v+width*k for v,k in zip(upper,response(lo,hi,True))]
        lo=hi
    degree = pot["short"]["base"]["arc_budget"]()[1].upper()/2
    mult = 32*arb.pi()**2*r**4
    tail = 608*degree**3/lo
    intervals = [v.lower().union((w+n*tail).upper())/mult for v,w,n in zip(lower,upper,(1,1,2,2))]
    a,b,plus,minus = intervals
    x = ((plus-minus)/4).intersection((plus-a-b)/2).intersection((a+b-minus)/2)
    x *= pot["m"]["axis"]["cos_axis"]
    det = a*b-x*x
    assert det>0
    factor = 1/(arb.pi()*det.sqrt())
    for name,val in zip(("H11","H22","H12","low-energy density slope"),(a,b,x,factor)):
        print(name,val.lower().str(12),val.upper().str(12))

if __name__=="__main__":
    check()
