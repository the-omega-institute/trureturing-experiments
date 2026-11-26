"""Mass-ratio floor for the Gaussian cross-ring covariance."""
from pathlib import Path
import runpy
from functools import lru_cache
from flint import acb,acb_mat,arb

work=Path(__file__).resolve().parent
shared=runpy.run_path(str(work/"ch39-schur-ring-hierarchy.py"))
s=shared["s"]
shell=s["shell"]
block=s["block"]
omega=shared["omega"]
wall=(6,7,8,11,12,18)
inside=tuple(i for i in range(19) if i not in wall)
modes=[(0,0),(-1,0),(0,-1)]
harmonics=[(1,0),(0,1),(1,-1)]
matrix_u=shared["shared"]["matrix_positive_parameter"]


def poly(m):
    mixed=block(m,wall,inside)
    pieces=[mixed,acb_mat(6,13),acb_mat(6,13)]
    for i,j,d,_ in shell["edges"]:
        if d in (724,804):
            assert i in inside and j in wall[:3]
            row,col=wall.index(j),inside.index(i)
            pieces[1+(d==804)][row,col]=mixed[row,col]
            mixed[row,col]=0
    inv=block(m,inside,inside).inv()
    data={(0,0):block(m,wall,wall)}
    for p,a in zip(modes,pieces):
        for q,b in zip(modes,pieces):
            k=p[0]-q[0],p[1]-q[1]
            data[k]=data.get(k,acb_mat(6,6))-a*inv*b.conjugate().transpose()
    return data


def symbols(m):
    data=poly(m)
    den,tau=[],[]
    for l in range(3):
        row=lambda z,i,j:sum(z[i,j+h]*omega**(h*l) for h in range(3))
        aa={k:row(z,0,0) for k,z in data.items()}
        bb=row(data[0,0],3,3).real
        assert bb>0
        t=[row(data[k],0,3) for k in modes]
        d={k:v*bb for k,v in aa.items()}
        for p,v in zip(modes,t):
            for q,w in zip(modes,t):
                key=p[0]-q[0],p[1]-q[1]
                d[key]-=v*w.conjugate()
        den.append(d)
        tau.append(t)
    return den,tau


def real_terms(d):
    return (d[0,0].real,[(2*d[key].real,-2*d[key].imag) for key in harmonics])


def positive(d,parts):
    center,terms=real_terms(d)
    total=sum((a**2+b**2).abs_upper().sqrt() for a,b in terms)
    if center>total:
        return True
    # Minimize over q exactly and enclose the remaining one-dimensional minimum.
    for j in range(parts):
        unit=unit_interval(j,parts)
        first=center+2*(d[1,0]*unit).real
        second=2*abs(d[0,1]+d[1,-1].conjugate()*unit.conjugate())
        if not first>second:
            return False
    return True


@lru_cache(maxsize=512)
def unit_interval(j,parts):
    angle=arb.pi()*arb(2*j+1,1)/parts
    return acb(angle.cos(),angle.sin())


def box(low,high,det_ratio,n):
    try:
        ds,rows=symbols(matrix_u(arb(low).union(high)))
    except (ZeroDivisionError,AssertionError):
        return False
    for i,d in enumerate(ds):
        if not positive(d,n):
            return False
        for j,other in enumerate(ds):
            if i!=j and not positive({k:v-det_ratio*other[k] for k,v in d.items()},n):
                return False
    for l,c in enumerate(rows):
        r,t=rows[(l+1)%3],rows[(l+2)%3]
        minors=[r[(j+1)%3]*t[(j+2)%3]-r[(j+2)%3]*t[(j+1)%3] for j in range(3)]
        ceiling=arb(".0086")*sum(map(abs,c))/det_ratio
        if not any(abs(abs(minors[(j+1)%3])-abs(minors[(j+2)%3]))
                   >ceiling*(abs(r[j])+abs(t[j])) for j in range(3)):
            return False
    return True


def check(low,high):
    ratio=(arb(".413")+arb(".045")*low).min(arb(".68"))
    if box(low,high,ratio,64):
        return 1
    assert high-low>arb(".00001"),(low,high,ratio)
    mid=(arb(low)+high)/2
    return check(low,mid)+check(mid,high)


if __name__=="__main__":
    count=check(arb(0),arb(6400))
    u=arb(6400)
    eta=arb(13)/(u+4-12)
    for i in range(19):
        degree=sum(i in (a,b) for a,b,_,_ in shell["edges"])
        assert shell["V"][i]+degree >= (4 if i in inside else 0)
    m=shell["matrix"](u,{})
    phase=(-m[6,12]).arg()
    ratio=abs(m[6,11])/abs(m[6,12])
    transfer=(ratio-1+arb(3).sqrt()*phase.sin())/(ratio+2)
    spread=(2+3*eta)/u
    det_ratio=((1-spread)**2-(3*(1+eta)/u)**2)/(1+spread)**2
    assert det_ratio>arb(".998") and det_ratio*(transfer-eta)/(1+eta)>arb(".020")
    mu=arb(".165")/(1-arb(".165"))
    normalized=arb(".0086")*(1-2*mu)/(1+2*mu)
    ceiling=arb(".00316")/arb(".8573")
    assert normalized>ceiling
    print("cross-covariance spectral intervals, ratio after point normalization, charm/top maximum",
          count,normalized.str(12),ceiling.str(12))
