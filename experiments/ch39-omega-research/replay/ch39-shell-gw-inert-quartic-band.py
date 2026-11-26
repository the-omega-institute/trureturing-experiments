"""A positive scalar potential through the first two extra thresholds in the 43rd-zero window."""
from itertools import combinations
from pathlib import Path
from runpy import run_path

from flint import arb, fmpq, fmpq_mat

folder = Path(__file__).resolve().parent
matched = run_path(str(folder/"ch39-shell-gw-smdr-matched-crossing.py"))
model,flow,q,loop = (matched[k] for k in ("model","flow","q","loop"))


def scalar_terms(diag,a,b):
    L = [24*l*l for l in diag]
    A,B = {},{}
    for (i,j),x in a.items():
        y = b[i,j]
        add = 2*x*x+2*x*y+y*y
        L[i] += add
        L[j] += add
        A[i,j] = 4*(diag[i]+diag[j])*(3*x+y)+4*x*x+2*y*y
        B[i,j] = 4*(diag[i]+diag[j])*y+8*x*y+4*y*y
        for k in range(len(diag)):
            if k in (i,j):
                continue
            ik,jk = sorted((i,k)),sorted((j,k))
            v,w = a[tuple(ik)],a[tuple(jk)]
            c,d = b[tuple(ik)],b[tuple(jk)]
            A[i,j] += 4*v*w+2*v*d+2*c*w
            B[i,j] += 2*c*d
    return L,A,B


def hessian_check():
    # Derive each spectator coefficient as a polynomial identity by polarization.
    # With three real four-vectors, x_i = |phi_i|² and the real/imaginary parts
    # of phi_i^dag phi_j are quadratic forms on the canonically normalized fields.
    n = 3
    pairs = list(combinations(range(n),2))
    layout = [("l",i,i) for i in range(n)]+[(t,i,j) for t in ("a","b") for i,j in pairs]
    identity = [[int(i==j) for j in range(4)] for i in range(4)]
    J = [[0,1,0,0],[-1,0,0,0],[0,0,0,1],[0,0,-1,0]]
    def form(i,j,block):
        M = fmpq_mat(12,12)
        for k in range(4):
            for l in range(4):
                M[4*i+k,4*j+l] = fmpq(block[k][l],2 if i==j else 4)
                if i!=j:
                    M[4*j+l,4*i+k] = M[4*i+k,4*j+l]
        return M
    norms = [form(i,i,identity) for i in range(n)]
    bi = [(t,i,j,[(norms[i],norms[j])] if t!="b" else [(S,S) for S in (form(i,j,identity),form(i,j,J))])
          for t,i,j in layout]
    def hess(x,S,T):
        sx,tx = S*x,T*x
        return 2*(x.transpose()*sx)[0,0]*T+2*(x.transpose()*tx)[0,0]*S+4*(sx*tx.transpose()+tx*sx.transpose())
    def flatten(H):
        return [z for z in H.entries()]
    def field(indices):
        v=fmpq_mat(12,1)
        for i,k in indices:
            v[4*i+k,0]=1
        return v
    samples = [[(i,0)] for i in range(n)]+[[(i,0),(j,k)] for k in (2,0) for i,j in pairs]
    hs = [[flatten(sum((hess(field(sample),S,T) for S,T in terms),fmpq_mat(12,12)))
           for _,_,_,terms in bi] for sample in samples]
    def verify(coeffs):
        L = coeffs[:3]
        A = dict(zip(pairs,coeffs[3:6]))
        B = dict(zip(pairs,coeffs[6:]))
        l,a,b = scalar_terms(L,A,B)
        expected = list(l)+[l[i]+l[j]+a[i,j] for i,j in pairs]+[l[i]+l[j]+a[i,j]+b[i,j] for i,j in pairs]
        for tensors,val in zip(hs,expected):
            combined = [sum(c*z for c,z in zip(coeffs,col)) for col in zip(*tensors)]
            assert 2*sum(v*v for v in combined)==val
    for i in range(9):
        verify([int(j==i) for j in range(9)])
        for k in range(i):
            verify([int(j in (i,k)) for j in range(9)])


def advance(k,duration,g,quart,y_light):
    n = k+1
    pairs = list(combinations(range(n),2))
    n_pairs=len(pairs)
    dt = duration/64
    slopes = [b+k*d for b,d in zip(flow["beta"],model["coeff"]["db"])]
    Bmat = [[a+k*e for a,e in zip(row,extra)] for row,extra in zip(model["coeff"]["base"],model["coeff"]["extra"])]
    for _ in range(64):
        assert 10*dt*y_light<1
        ymax = y_light/(1-10*dt*y_light)
        gs = [v+arb(0,".009") for v in g]
        assert min(gs)>0
        dg = []
        for i,b in enumerate(slopes):
            D = model["coeff"]["yukawa"][i]
            damp = q(0).union(max(D)*ymax)
            dg.append(2*gs[i]**2*(b+(sum(x*y for x,y in zip(Bmat[i],gs))-damp)/loop))
        assert all(abs(dt*v)<q(".0089") for v in dg)
        W,Y = gs[1],q(3,5)*gs[0]
        # Only the light doublet couples to fermions. T <= 3*ymax and
        # Tr(3U²+3D²+E²) <= 3*ymax².
        trace,trace4 = q(0).union(3*ymax),q(0).union(3*ymax**2)
        box = [v+arb(0,".009") for v in quart]
        L = box[:n]
        a = dict(zip(pairs,box[n:n+n_pairs]))
        b = dict(zip(pairs,box[n+n_pairs:]))
        dl,da,db = scalar_terms(L,a,b)
        common = 3*(3*W+Y)
        for i in range(n):
            dl[i] += q(3,8)*(2*W**2+(W+Y)**2)-common*L[i]
        dl[0] += 4*trace*L[0]-2*trace4
        for i,j in pairs:
            shift = (2*trace if i==0 else 0)-common
            da[i,j] += q(3,4)*(3*W**2+Y**2-2*W*Y)+shift*a[i,j]
            db[i,j] += 3*W*Y+shift*b[i,j]
        beta = [*dl,*da.values(),*db.values()]
        assert all(abs(dt*v)<q(".0089") for v in beta)
        tube = [v+q(0).union(dt)*d for v,d in zip(quart,beta)]
        assert min(tube[:n])>0
        assert all(tube[n+i]>0 and tube[n+i]+tube[n+n_pairs+i]>0 for i in range(n_pairs))
        quart = [v+dt*d for v,d in zip(quart,beta)]
        g = [v+dt*d for v,d in zip(g,dg)]
        y_light = ymax
    return g,quart,y_light


def check():
    hessian_check()
    _,m = flow["masses"](43)
    eps = q("12e-8")/((m[0]**2-125**2)/q("1.2e12")**2)
    b = ((1+eps)/(1-eps)).sqrt()
    spread = (1/b).union(b)
    s0,s1 = q(".36"),q(".391")
    masses = [v*spread for v in m[:3]]
    first = model["shifted"](masses[0],s0.union(s1))
    ys,_ = matched["get_boundary"](first.lower(),first.upper())
    L = ys[6]
    assert L>q(".00091")
    g = [q(5,3)*ys[2]**2,ys[1]**2,ys[0]**2]
    y_light = sum(t*t for i,t in enumerate(ys) if i not in (0,1,2,6))
    # Seed both extra self-couplings and their mutual portals. A strictly
    # positive mixed exchange coupling also makes the positivity test simple.
    quart = [L,q(".10"),q(".60"),q(".01")]
    for k in (1,2):
        def distance(s):
            return (model["shifted"](masses[k],s)/model["shifted"](masses[k-1],s)).log()/loop
        duration = distance(s0).union(distance(s1))
        assert duration>0
        if k==2:
            l0,l1,a01,b01 = quart
            quart = [l0,l1,q(".10"),a01,q(".60"),q(".10"),b01,q(".01"),q(".01")]
        g,quart,y_light = advance(k,duration,g,quart,y_light)
    assert quart[0]>q(".01") and max(quart)<1
    print("quartic at the third threshold",float(quart[0].lower()),float(quart[0].upper()),
          "largest component upper",max(float(v.upper()) for v in quart))


if __name__=="__main__":
    check()
