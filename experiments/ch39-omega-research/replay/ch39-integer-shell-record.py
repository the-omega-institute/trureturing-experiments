"""Replay the shell, arithmetic Gram and Schur mass with Python integer arithmetic."""
import json
from fractions import Fraction as Q
from itertools import combinations, product
from pathlib import Path


def turn(chart, den):
    s,n=chart
    assert s in (-1,1)
    return Q(s*(den**2-n*n),den**2+n*n),Q(s*2*n*den,den**2+n*n)


def times(a,b):
    return a[0]*b[0]-a[1]*b[1],a[0]*b[1]+a[1]*b[0]


def main():
    data=json.loads(Path(__file__).with_name("ch39-integer-shell-record.json").read_text())
    den=2**data["bits"]
    rows=data["rows"]
    assert all(r[2]>0 for r in rows) and sum(r[2] for r in rows)==den
    phase=[Q(0)]*4
    f=[Q(0)]*4
    for _,height,num,first,second in rows:
        g,w=Q(height,den),Q(num,den)
        a,b=turn(first,den),turn(second,den)
        assert all(x*x+y*y==1 for x,y in [a,b])
        values=[a[0],b[0],times((a[0],-a[1]),b)[0],times(a,b)[0]]
        phase=[v+w*y for v,y in zip(phase,values)]
        f=[v+w*g*g/(g*g+j*j) for j,v in enumerate(f)]
    assert f[0]==1
    e=max(map(abs,phase))
    assert e<Q(1,2**204)
    # For the nine shell states coefficient l1 norm is <=2, so the
    # normalized inner products differ from their exact targets by <=8e/(1-4e).
    points=[v for v in product(range(5),repeat=3) if sum(i*i for i in v)==17]
    target=lambda u,v:Q(17+sum(a*b for a,b in zip(u,v)),34)
    assert len(points)==9 and max(target(u,v) for u,v in combinations(points,2))==Q(33,34)
    assert target((0,1,4),(2,3,2))==Q(14,17) and 8*e/(1-4*e)<Q(1,2**200)
    def inner(i,j):
        return f[i] if i==j else (i*f[i]+j*f[j])/(i+j)
    for j in range(1,4):
        for i in range(j):
            raw=inner(i,j)
            h2=raw*raw/(f[i]*f[j])
            k=data["gram_bins"][j][i]
            assert raw>0 and Q(k,2**96)**2<h2<Q(k+1,2**96)**2
    a,b,c=inner(0,1),inner(1,2),inner(0,2)
    minor=f[0]*f[2]-c*c
    det=f[0]*f[1]*f[2]+2*a*b*c-f[0]*b*b-f[1]*c*c-f[2]*a*a
    assert det>0 and Q("528.89869680")**2<f[1]*minor/det<Q("528.89869682")**2
    print("Integer replay: 9 shell points, arc targets 33/34 and 14/17 within 2^-200.")
    print("Six 96-bit Gram bins and 528.89869680 < M_F/M_S < 528.89869682.")


if __name__=="__main__":
    main()
