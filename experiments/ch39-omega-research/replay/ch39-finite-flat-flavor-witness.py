"""Finite rational flat-site rotors giving six quark eigenvalues and CKM data in the light ground."""
import json
from itertools import permutations
from pathlib import Path

from flint import arb, acb, acb_mat, ctx

folder = Path(__file__).resolve().parent
q = arb


def cs(t):
    return (1-t*t)/(1+t*t),2*t/(1+t*t)


def build():
    raw = json.loads((folder/"ch39-finite-flat-rotors.json").read_text())
    assert len(raw["rotors"])==24 and raw["g"]==5
    pts = list(permutations((0,1,3,4)))
    a,b = pts[:2]
    alpha = (q(5).sqrt()-1)/2
    labels = dict(zip(a,(1,alpha,-alpha,-1)))
    up,down = acb_mat(3,3),acb_mat(3,3)
    for p,(sign,integer) in zip(pts,raw["rotors"]):
        assert sign in (1,-1) and abs(integer)<=raw["den"]
        r,i = cs(q(integer)/raw["den"])
        z = sign*acb(r,i)
        v = [labels[t] for t in p]
        col = acb_mat([[v[0]+v[j]] for j in (1,2,3)])
        site = col*col.transpose()/(8*(2+2*alpha**2))
        chi = 1 if sum(t*(x-y) for t,x,y in zip(p,a,b))>0 else -1
        h = 5*z/q(24).sqrt()
        up += h*site
        down += chi*h.conjugate()*site
    return up,down


def projectors(m,values):
    H = m*m.conjugate().transpose()
    poly = H.charpoly()
    roots = []
    for y in values:
        lo,hi = (y-q("1e-13"))**2,(y+q("1e-13"))**2
        assert (poly(lo).real*poly(hi).real) < 0
        z = lo.union(hi)
        for _ in range(5):
            x = z.mid()
            z = z.intersection(x-poly(x).real/poly.derivative()(z).real)
        roots.append(z)
    assert all(x<y for x,y in zip(roots,roots[1:]))
    I = acb_mat([[int(i==j) for j in range(3)] for i in range(3)])
    out = []
    for j,z in enumerate(roots):
        P = I
        for k,t in enumerate(roots):
            if k!=j:
                P = P*(H-t*I)/(z-t)
        out.append(P)
    return out


def main():
    ctx.prec = 320
    up,down = build()
    P = projectors(up,[q(v) for v in (".00000615",".00311",".8616")])
    Q = projectors(down,[q(v) for v in (".0000135",".000268",".01401")])
    a,b = cs(q(57)/500)
    c,d = cs(q(107)/5000)
    e,f = cs(q(377)/200000)
    re,im = cs(q(16)/25)
    phase = acb(re,im)
    V = acb_mat([[1,0,0],[0,c,d],[0,-d,c]])*acb_mat(
        [[e,0,f*phase.conjugate()],[0,1,0],[-f*phase,0,e]])*acb_mat([[a,b,0],[-b,a,0],[0,0,1]])
    tr = lambda m:sum(m[i,i] for i in range(3))
    err = max(abs(tr(p*s).real-abs(V[i,j])**2).upper() for i,p in enumerate(P) for j,s in enumerate(Q))
    J = tr(P[0]*Q[0]*P[1]*Q[1]).imag
    assert err<q("1e-9") and abs(J-a*c*e**2*b*d*f*im)<q("1e-11")
    print("Mass singular errors < 1e-13; largest CKM squared-modulus error",err.str(8),
          "and signed J",J.str(12))


if __name__=="__main__":
    main()
