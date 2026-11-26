"""Mass hierarchy for transfer directly from the inner ring to the outer ring."""
from pathlib import Path
import runpy
from flint import acb,acb_mat,arb

shared=runpy.run_path(str(Path(__file__).with_name("ch39-schur-maxwell-mixing.py")))
s=shared["s"]
shell=s["shell"]
interior=s["inside"]
block=s["block"]
omega=acb(-arb(1)/2,arb(3).sqrt()/2)
threshold=arb(".013")


def fourier_rows(m):
    cols=(11,12,18)
    inv_rhs=block(m,interior,interior).solve(block(m,interior,cols))
    pieces=[acb_mat(1,len(interior)) for _ in range(3)]
    for i,site in enumerate(interior):
        edge=s["geometry"]["edges"].get(tuple(sorted((6,site))))
        if edge is None:
            continue
        distance=edge[0]
        k=1 if distance==724 else 2 if distance==804 else 0
        pieces[k][0,i]=m[6,site]
    bonds=[p*inv_rhs for p in pieces]
    bonds[0]-=block(m,(6,),cols)
    return [[sum(bond[0,j]*omega**(l*j) for j in range(3)) for bond in bonds] for l in range(3)]


def verified(low,high):
    try:
        coeff=fourier_rows(shared["matrix_positive_parameter"](arb(low).union(high)))
    except ZeroDivisionError:
        return False
    for l in range(3):
        r,t=coeff[(l+1)%3],coeff[(l+2)%3]
        upper=sum(abs(x) for x in coeff[l])
        minors=[r[(j+1)%3]*t[(j+2)%3]-r[(j+2)%3]*t[(j+1)%3] for j in range(3)]
        if not any(abs(abs(minors[(j+1)%3])-abs(minors[(j+2)%3]))
                   > threshold*upper*(abs(r[j])+abs(t[j])) for j in range(3)):
            return False
    return True


def verify_range(low,high):
    if verified(low,high):
        return 1
    assert high-low>arb(".00001")
    mid=(arb(low)+high)/2
    return verify_range(low,mid)+verify_range(mid,high)


if __name__=="__main__":
    count=verify_range(0,6400)
    m=shell["matrix"](arb(6400),{})
    d=-m[6,12]
    r=-m[6,11]
    assert r.real>0 and r.imag==0
    phi=d.arg()
    assert 0<phi<arb(".03") and abs(r)>abs(d)
    # The Euler product makes f(long)/f(short) decreasing.
    ratio=abs(r)/abs(d)
    tail=(ratio-1+arb(3).sqrt()*phi.sin())/(ratio+2)
    eta=arb(12)/(6400+4-11)
    assert (tail-eta)/(1+eta)>arb(".020")
    mu=arb(".165")/(1-arb(".165"))
    normalized=threshold*(1-2*mu)/(1+2*mu)
    up=arb(".00316")/arb(".8573")
    print("spectral parameter intervals, bare ratio floor",count,threshold.str(5))
    print("normalized transfer floor and input charm/top ceiling",normalized.str(12),up.str(12))
    assert normalized>arb(".00563") and normalized>up
