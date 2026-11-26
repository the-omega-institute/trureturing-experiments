"""Fit the same-scale quark data inside the first five singlet/adjoint scalar modes."""
from pathlib import Path
from runpy import run_path
import numpy as np
from flint import acb,acb_mat,arb


folder=Path(__file__).resolve().parent
m=run_path(str(folder/"ch39-shell-gw-higgs-low-band.py"))
yuk=m["yuk"];adj=m["adj"];hs=m["hs"]
m["ctx"].prec=400
up=[
    '-0.0034247194749082678 0.017531674736681786',
    '0.053414069551493068 0.0045584804774757492',
    '-0.00081522633671574443 0.0033188665674263873',
    '0.24534189920375035 0.097509303144303769',
    '0.2136833549045305 -0.73558341310402164',
    '0.038709434465304768 0.018355460162871966',
    '-0.025722841130424376 0.086247522567274471',
    '0.2641162594221827 0.052365402960649808',
    '-0.0063882915225605642 0.015472205656220943']
down=[
    '5.4484780007049607e-05 0.00032967992871629401',
    '0.00089310224776674278 -0.00012037128874981155',
    '8.2952768685543242e-06 6.6013725557971051e-05',
    '0.0044313042917422454 0.00068434440378722049',
    '0.0024000384617240719 -0.012070043784355812',
    '0.00071555436755990122 0.00015762415231806016',
    '-9.8729886514854604e-05 0.0015238820634413922',
    '0.0045647530904753706 -0.00012634987953519621',
    '0.00011401204550494705 0.00037177058643298588']


def refined_frame(k,n):
    rows=[[k[i,j] for j in range(k.ncols())] for i in range(k.nrows())]
    count=m["model"]["model"]["count"]
    cut=arb("1.087")
    assert count(rows,cut)==n and count(rows,arb("1.0871"))==n
    vals,eigen=np.linalg.eigh(np.array([[complex(z) for z in row] for row in rows]))
    result=[]
    for i in range(n):
        v=acb_mat([[acb(float(z.real),float(z.imag))] for z in eigen[:,i]])
        e=arb(float(vals[i]))
        # Correct the approximate eigenvector while keeping the solve away from singularity.
        shift=e+acb(0,arb("1e-11"))
        matrix=k-shift*yuk["identity"](k.nrows())
        for _ in range(5):
            v=matrix.solve(v)*(-acb(0,arb("1e-11")))
        e=(hs(v,k*v)/hs(v,v)).real
        err=hs(k*v-e*v,k*v-e*v).real.sqrt()/(cut-e)
        assert err<arb("1e-23"), (n,i,err,hs(v,v))
        margin=arb(0,err.upper())
        result.append(acb_mat([[v[j,0]+acb(margin,margin)] for j in range(k.nrows())]))
    return result


def matrices():
    sing,fam,phases,g=m["transported_spectra"]()
    L,R=yuk["zero_basis"](1),yuk["zero_basis"](-1)
    sources=[(x,[yuk["identity"](3)/arb(3).sqrt()]) for x in refined_frame(sing,1)]
    sources.extend((x,g) for x in refined_frame(fam,4))
    out=[]
    for v,G in sources:
        y=acb_mat(3,3)
        for i in range(6):
            H=sum((v[i*len(G)+b,0]*t for b,t in enumerate(G)),acb_mat(3,3))
            a=acb_mat([[L[3*i+r,c] for c in range(3)] for r in range(3)])
            b=acb_mat([[R[3*i+r,c] for c in range(3)] for r in range(3)])
            y+=phases[i]*a.transpose()*H.transpose()*b
        out.append(y)
    return sources,out


def fit(target,source,ys):
    t=acb_mat([[acb(*(arb(z) for z in target[3*i+j].split())) for j in range(3)] for i in range(3)])
    gram=acb_mat([[hs(x,y) for y in ys] for x in ys])
    coef=gram.solve(acb_mat([[hs(y,t)] for y in ys]))
    y=sum((coef[i,0]*z for i,z in enumerate(ys)),acb_mat(3,3))
    err=hs(y-t,y-t).real.sqrt()
    assert err<arb("1e-10")
    norm=acb_mat([[hs(a,b) if len(G)==len(H) else 0 for b,H in source] for a,G in source])
    return y,(adj(coef)*norm*coef)[0,0].real.sqrt(),err


def spec(y):
    h=y*adj(y);vs=[z.real for z in h.charpoly().roots(tol=1e-19)]
    vs.sort(key=lambda a:float(a.mid()))
    ps=[]
    for i,e in enumerate(vs):
        p=yuk["identity"](3)
        for j,f in enumerate(vs):
            if j!=i:p=p*(h-f*yuk["identity"](3))/(e-f)
        ps.append(p)
    return [e.sqrt() for e in vs],ps


if __name__=="__main__":
    op,ys=matrices()
    upy,upnorm,eu=fit(up,op,ys)
    downy,dnorm,ed=fit(down,op,ys)
    u,p=spec(upy);d,q=spec(downy)
    targets=[arb(a) for a in (".00000615",".00311",".8616",".0000135",".000268",".01401")]
    assert all(abs(v/ref-1)<arb("1e-6") for v,ref in zip(u+d,targets))
    tr=m["trace"];vs=[[tr(a*b).real for b in q] for a in p]
    z=vs[0][2];s13=z.sqrt();c13=(1-z).sqrt()
    s12=vs[0][1].sqrt()/c13;s23=vs[1][2].sqrt()/c13
    quartet=tr(p[0]*q[1]*p[1]*q[2])
    delta=(quartet+z*s12**2*s23**2*c13**2).arg()
    angles=[s12.asin(),s23.asin(),s13.asin(),delta]
    assert all(abs(v-arb(t))<arb("1e-6") for v,t in zip(angles,(".22704",".04275",".00377","1.139")))
    print("five modes below 1.087, scalar norms",upnorm.str(10),dnorm.str(10),
          "fitted J",quartet.imag.str(10))
