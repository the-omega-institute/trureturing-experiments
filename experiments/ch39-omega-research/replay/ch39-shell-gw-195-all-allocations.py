"""Classify both sector bandwidths after the first Wilson-conditioned zero switch."""
from pathlib import Path
from runpy import run_path
import numpy as np
from flint import acb, acb_mat, arb

m=run_path(str(Path(__file__).with_name("ch39-shell-gw-195-flavor.py")))
four=m["bands"]["four"]
adj,hs=m["adj"],m["hs"]
a=arb


def texture(B):
    diagonal=acb_mat([[b[i,i] for b in B[:3]] for i in range(3)])
    highdiag=acb_mat([[b[i,i] for b in B[3:5]] for i in range(3)])
    low={}
    high={}
    for i in range(3):
        for j in range(3):
            low[i,j]=acb_mat([[b[i,j] for b in B[:3]]])*diagonal.inv()
            high[i,j]=acb_mat([[b[i,j] for b in B[3:5]]])-low[i,j]*highdiag
    return low,high


def adapt(v,p):
    trial=np.eye(3,dtype=complex)
    trial[:,p]=v
    order=[p]+[i for i in range(3) if i!=p]
    Q,_=np.linalg.qr(trial[:,order])
    return four["exact_frame"](Q[:,np.argsort(order)])


def certify_textures(B):
    L=four["unitary_trial"](np.linalg.eig(m["array"](B[1])@np.linalg.inv(m["array"](B[0])))[1])
    R=four["unitary_trial"](m["array"](B[0]).conj().T@L)
    L,R=four["exact_frame"](L),four["exact_frame"](R)
    m["four_exclusion"](B[:4], (L,R))
    B=[adj(L)*b*R for b in B[:5]]
    old,_=texture(B)
    def stack(dic,positions):
        return acb_mat([[dic[i,j][0,k] for k in range(dic[i,j].ncols())] for i,j in positions])
    for p in range(3):
        i,j=[k for k in range(3) if k!=p]
        first=acb_mat([[old[r,c][0,p] for c in range(3)] for r in range(3)])
        left,_,right=np.linalg.svd(m["array"](first))
        ul,ur=adapt(left[:,0],p),adapt(right.conj().T[:,0],p)
        for u in (ul,ur):
            assert (1-abs(u[p,p])**2).sqrt()<a(".0011")
        low,high=texture([adj(ul)*b*ur for b in B])
        for r,c in ((i,j),(j,i),(i,p),(j,p),(p,i),(p,j)):
            assert all(abs(low[r,c][0,k]) < a(".00135") for k in range(3))
            if p in (r,c):
                assert abs(low[r,c][0,p])<a("1e-8")
        for pos in ([(i,p),(j,p)],[(p,i),(p,j)]):
            T=stack(high,pos)*stack(high,[(i,j),(j,i)]).inv()
            test=a("3.46")**2*m["eye"](2)-adj(T)*T
            assert test[0,0].real>0 and test.det().real>0
            f=stack(low,pos)-T*stack(low,[(i,j),(j,i)])
            for k in range(3):
                limit=".00130" if k==p else ".00366"
                assert sum(abs(f[h,k])**2 for h in range(2)).sqrt()<a(limit)
        for iteration in range(2):
            T=stack(high,[(i,p),(j,p)])*stack(high,[(i,j),(j,i)]).inv()
            t=T[0,1]/T[1,1]
            if abs(t)<1:break
            i,j=j,i
        f=stack(low,[(i,p),(j,p)])-T*stack(low,[(i,j),(j,i)])
        assert a(".233")<abs(t)<a(".234") and abs(T[0,0]-t*T[1,0])<a(".234")
        assert abs(f[0,p]-t*f[1,p])<a(".000007")
        assert max(abs(f[0,k]-t*f[1,k]) for k in (i,j))<a(".00124")
        assert abs(low[i,j][0,p])<a(".000029") and abs(low[j,i][0,p])<a(".000388")


def transfers(x,r,first):
    # If either right complementary coordinate is unknown, each other diagonal
    # entry is bounded initially by x+r, later by x*y+r.
    D=a(".00130")+2*a(".00366")*(x+r)
    y=a(1)
    for _ in range(4):
        fixed=a("3.46")*(r**2+first**2).sqrt()+D
        y=(fixed+r*x)/((1-x*x).sqrt()-a("3.46")*x)
        assert y<1
        xnew=(a("3.46")*x*y+fixed+r*y)/(1-y*y).sqrt()
        x=min(x,xnew)
        D=a(".00130")+2*a(".00366")*(x*y+r)
    return x,y


def test_cases():
    tol=a("1e-6")
    t12,t23,t13=[a(v)+a(0,tol) for v in (".22704",".04275",".00377")]
    vub,vcb=t13.sin(),t23.sin()*t13.cos()
    mix=(1-(t23.cos()*t13.cos())**2).sqrt()
    ru,rd,fu,fd=a(".00362"),a(".01914"),a(".0000072"),a(".00097")
    mult=(1+tol)/(1-tol)
    assert a(".00311")/a(".8616")*mult<ru and a(".000268")/a(".01401")*mult<rd
    assert a(".00000615")/a(".8616")*mult<fu and a(".0000135")/a(".01401")*mult<fd
    assert a(".00311")/a(".8616")/mult>a(".00360")
    assert a(".0428")<mix<a(".043")
    # Four down modes: compare the adapted axis to the original pivot.
    down=a(".02428")+a(".0011")
    xu,yu=transfers(mix+down,ru,fu)
    assert xu < mix-down
    # Three up modes: first obtain diagonal dominance after adapting the pivot.
    raw=a(2).sqrt()*a(".00131")/(1-ru-4*a(".00131"))
    x=raw+a(".0011")
    otherdiag=x*x+ru
    assert otherdiag<a(".00364")
    off=a("1e-8")+2*a(".00135")*otherdiag
    # Both cross vectors have norm <=sqrt(2)*off and the complement block
    # has norm <=x²+ru by the rank-one and remainder decomposition.
    top=a(2).sqrt()*off/(1-otherdiag)
    assert top<a(".000015")
    mu=a(".00360")
    ab=(a(".000029")+2*a(".00135")*otherdiag+top**2+fu)/mu
    ba=(a(".000388")+2*a(".00135")*otherdiag+top**2+fu)/mu
    large=((1-top**2)/2).sqrt()
    _,beta=transfers(mix+top,rd,fd)
    assert beta<a(".09")
    for light_i in (True,False):
        direct,reverse=(ab,ba) if light_i else (ba,ab)
        charm=direct/(1-top**2-(reverse/large)**2).sqrt()
        light=(top**2+charm**2).sqrt()
        small=vub+charm*vcb+top
        big=(1-light**2).sqrt()*vcb-light*vub-top
        si=t12.sin()*t13.cos()+charm+top
        row=((rd*si)**2+fd**2).sqrt() if light_i else rd
        bi=small if light_i else mix+top
        lhs=(a(".233")*big-small if light_i else big-a(".234")*small)*(1-beta**2).sqrt()
        error=a(".000007")+2*a(".00124")*((mix+top)*beta+rd)
        rhs=a(".234")*(bi*beta+row)+(row+a(".234")*rd)*beta+error
        assert lhs>rhs
        print("three-up five-down mismatch margin",(lhs-rhs).lower().str(8))
    print("five-up four-down axis discrepancy",(mix-down-xu).lower().str(8))


if __name__=="__main__":
    sg,fa,F=m["scalar_models"]()
    B=m["image_basis"](F,m["sources"](sg,fa,5))
    certify_textures(B)
    test_cases()
    from itertools import combinations
    pairs=list(combinations(range(3),2))
    def wedge(x,y):
        return acb_mat([[x[i,r]*y[j,s]-x[i,s]*y[j,r] for r,s in pairs] for i,j in pairs])
    b,c=B[:2]
    cols=[wedge(b,b),(wedge(b,c)+wedge(c,b))/a(2).sqrt(),wedge(c,c)]
    G=acb_mat([[hs(x,y) for y in cols] for x in cols])
    assert m["count"](m["bands"]["m"]["rows"](G),a(".28")**2)==0
    r=a(".20")
    assert (2*r**2+r**4)/(1+r**2)**2<a(".28")**2
