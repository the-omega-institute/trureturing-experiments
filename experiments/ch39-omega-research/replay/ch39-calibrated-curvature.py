"""Certified positive-measure boxes at fixed volume, with every squared edge integrated."""
from itertools import combinations
from pathlib import Path
from runpy import run_path

from flint import arb, arb_mat, fmpq

switches = run_path(str(Path(__file__).with_name("ch39-calibrated-length-switches.py")))
angle, gram = switches["angle"], switches["gram"]
pairs = switches["pairs"]
volume_per_cell = 1/arb(5).sqrt()
eta = fmpq(1,2**40)


def local(parity, pos):
    q = {(i,i):arb(0) for i in range(5)}
    for i,j in pairs:
        d = j-i
        center = fmpq(d*(5-d),5)
        if parity:
            center *= (1+fmpq(parity*(-1)**i,24))*(1+fmpq(parity*(-1)**j,24))
        elif (i,j) == (pos,pos+1):
            center += fmpq(1,3)
        q[i,j] = q[j,i] = arb(center,eta)
    g = gram(q,range(5))
    for k in range(1,5):
        assert arb_mat([[g[i,j] for j in range(k)] for i in range(k)]).det() > 0
    det = g.det()
    assert all(det > v**4/16 for v in q.values())
    d = arb_mat([[-1]*4]+[[int(i==j) for j in range(4)] for i in range(4)])
    normals = d*g.inv()*d.transpose()
    s = arb(0)
    for i,j in pairs:
        face = tuple(v for v in range(5) if v not in (i,j))
        area = gram(q,face).det().sqrt()/2
        cos = -normals[i,j]/(normals[i,i]*normals[j,j]).sqrt()
        theta = angle(cos.lower()).union(angle(cos.upper()))
        degree = 6 if j-i in (1,4) else 4
        s += area*(fmpq(1,degree)-theta)
    return det.sqrt()/24,s


def verify():
    for sign,orbits in [(1,[(p,0,12) for p in (-1,1)]),
                        (-1,[(0,j,6) for j in range(4)])]:
        vol,s = arb(0),arb(0)
        for parity,pos,count in orbits:
            v,action = local(parity,pos)
            vol += count*v
            s += count*action
        scale = (volume_per_cell/vol).sqrt()
        assert scale > fmpq(1,2)
        favorable = sign*s*scale
        assert favorable > fmpq(1,250)
        print(sign,"signed action per cell at the flat volume",favorable,"scale",scale)
    assert 3*(960*volume_per_cell/15).sqrt() < 32
    energy = 30484476*9*(fmpq(1,250)-fmpq(1,400))
    assert switches["units"]() > 9
    assert energy-46*15 > 410000 and energy/256-46*15 > 900
    assert 1/(400**2*60*50*volume_per_cell) > fmpq(1,2**30)+4*fmpq(1,2**31)


if __name__ == "__main__":
    verify()
