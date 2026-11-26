"""Orbit-symmetric positive strengths realize target masses in degree zero, at RH zero eight."""
from pathlib import Path
import runpy
from flint import arb,acb,acb_mat,arb_mat

work=Path(__file__).resolve().parent
src=runpy.run_path(str(work/"ch39-single-arc-source.py"))
s,J,W=src["shell"],src["J"],src["W"]
block=src["block"]
even={0,3,4}
omega=acb(-arb(1)/2,arb(3).sqrt()/2)
inner,outer=src["inner"],src["outer"]
ii,oo=[[W.index(i) for i in arm] for arm in (inner,outer)]


def model(p,q):
    m=s["matrix"](0,{724:acb(0,p).exp(),804:acb(0,q).exp()})
    orbits=[[] for _ in range(4)]
    for i,j,d,z in s["edges"]:
        if d not in (724,804):
            continue
        v=acb_mat(19,1)
        v[i,0]=1
        v[j,0]=m[i,j].conjugate()
        orbits[2*(i not in even)+(d==804)].append(v)
    return m,orbits


def symbols(t):
    return [sum(t[0,j]*omega**(-l*j) for j in range(3)) for l in range(3)]


def evaluate(params,weights,goals):
    base,orbits=params
    m=acb_mat(base)
    for w,orbit in zip(weights,orbits):
        for h in orbit:
            m+=w*h*h.conjugate().transpose()
    c=block(m,J,J)
    projection=block(m,W,J)*c.inv()
    k=block(m,W,W)-projection*block(m,J,W)
    eig=symbols(-block(k,ii,oo))
    def components(z):
        delta=[z[l+1]-goals[l]*z[0] for l in range(2)]
        return [w.real for w in delta]+[w.imag for w in delta]
    columns=[]
    for orbit in orbits:
        derivative=acb_mat(3,3)
        for h in orbit:
            v=block(h,W,(0,))-projection*block(h,J,(0,))
            derivative-=block(v,ii,(0,))*block(v,oo,(0,)).conjugate().transpose()
        columns.append(components(symbols(derivative)))
    return arb_mat([[w] for w in components(eig)]),arb_mat(columns).transpose(),k


def start(raw,params):
    x,y,a,b=map(arb,raw)
    m,_=params
    g=block(m,J,J).inv()[0,0].real
    even_total=x/((1-x)*g)
    odd_total=y/((1-y)*g)
    return [even_total*a,even_total*(1-a),odd_total*b,odd_total*(1-b)]


def certify(raw,goals,exact,angular):
    center=start(raw,exact)
    for _ in range(4):
        residual,jacobian,k=evaluate(exact,center,goals)
        change=jacobian.solve(residual)
        center=[(v-change[i,0]).mid() for i,v in enumerate(center)]
    radius=arb(".000001")
    box=[v+arb(0,radius) for v in center]
    residual,j0,_=evaluate(angular,center,goals)
    inv=evaluate(exact,center,goals)[1].inv()
    _,j,k=evaluate(angular,box,goals)
    step=inv*residual
    perturb=arb_mat([[int(i==l) for l in range(4)] for i in range(4)])-inv*j
    for i in range(4):
        assert abs(step[i,0])+radius*sum(abs(perturb[i,l]) for l in range(4))<radius
    assert min(box)>0
    assert abs(symbols(-block(k,ii,oo))[0])>1
    loops={}
    for f in src["all_faces"]():
        z=acb(1)
        for a,b in zip(f,f[1:]+f[:1]):
            z*=-k[W.index(a),W.index(b)]
        loops[f]=[z]
    assert src["degrees"](loops,0)==(0,[0])
    print("positive stiffnesses, even 724/804 then odd 724/804",*[w.str(9) for w in center])


if __name__=="__main__":
    height=acb.zeta_zero(8).imag
    p,q=[s[key]*height for key in ("a","b")]
    exact=model(p,q)
    noise=arb(0,arb("1e-11"))
    angular=model(p+noise,q+noise)
    certify((".53791044",".42195445",".50656433",".64188365"),(0,0),exact,angular)
    certify((".53403701",".40691411",".50449069",".64941623"),
            (arb(".00311")/arb(".8616"),arb("6.15e-6")/arb(".8616")),exact,angular)
    certify((".51398282",".34135819",".49580459",".68872749"),
            (arb(".000268")/arb(".01401"),arb("1.35e-5")/arb(".01401")),exact,angular)
