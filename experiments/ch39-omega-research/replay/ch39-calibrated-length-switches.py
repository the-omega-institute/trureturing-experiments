"""Finite Regge length switches at the common Newton/gauge calibration."""
from itertools import combinations
from math import comb
from pathlib import Path
from runpy import run_path

from flint import arb, arb_mat, fmpq

source = run_path(str(Path(__file__).with_name("ch39-regge-fixed-volume.py")))
angle = source["angle_turns"]
vertices = tuple(range(5))
pairs = tuple(combinations(vertices,2))


def gram(q, face):
    start,*ends = face
    return arb_mat([[(q[start,i]+q[start,j]-q[i,j])/2 for i in ends] for j in ends])


def one_position(pos, shift):
    distances = {(i,i): fmpq(0) for i in vertices}
    for i,j in pairs:
        d = j-i
        q = fmpq(d*(5-d),5) + (shift if (i,j) == (pos,pos+1) else 0)
        distances[i,j] = distances[j,i] = q
    g = gram(distances,vertices)
    assert g.det().overlaps((1+shift-3*shift**2/4)/5)
    for k in range(1,5):
        assert arb_mat([[g[i,j] for j in range(k)] for i in range(k)]).det() > 0
    assert abs(shift) < fmpq(1,16) and g.det()/fmpq(6,5)**4 > fmpq(1,16)
    d = arb_mat([[-1]*4]+[[int(i==j) for j in range(4)] for i in range(4)])
    normals = d*g.inv()*d.transpose()
    action = arb(0)
    for i,j in pairs:
        face = tuple(v for v in vertices if v not in (i,j))
        area = arb(gram(distances,face).det()).sqrt()/2
        cosine = -arb(normals[i,j])/arb(normals[i,i]*normals[j,j]).sqrt()
        degree = 6 if j-i in (1,4) else 4
        action += area*(fmpq(1,degree)-angle(cosine))
    return action


def units():
    # Inscribed and circumscribed perimeters of the regular 1024-gon.
    cos = source["cosine_tick"](source["GRID"]//2048)
    sin = (1-cos*cos).sqrt()
    turn = arb(fmpq(6283,1000),fmpq(1,1000))
    assert 2048*sin > turn.lower() and 2048*sin/cos < turn.upper()
    # ln(2) = 2 atanh(1/3), with a geometric tail after the seventh term.
    lo = sum(fmpq(2,(2*j+1)*3**(2*j+1)) for j in range(7))
    hi = lo + fmpq(2,15*3**15)*fmpq(9,8)
    return turn/arb((lo+hi)/2,(hi-lo)/2)


def verify():
    # The original calibrated value and both adjacent integer choices lie in this hull.
    coupling = arb(fmpq(60968953,2),fmpq(1,2))
    to_bits = coupling*units()
    for n,lo,hi in [(6,7117,7125),(10,28,29),(11,7,8),(12,1,2),(13,0,1)]:
        shift = fmpq(1,2**n)
        s = 6*sum(one_position(i,shift) for i in range(4))
        bits = -to_bits*s
        assert lo < bits and bits < hi
        print("shift 2^-", n, "action in turns", s, "suppression bits", bits)
    costs = []
    for sign in (-1,1):
        h = fmpq(sign,128)
        x = 2*h+h*h
        shift = 2*x/(1+arb(1-3*x).sqrt())
        assert ((1+shift-3*shift**2/4).sqrt()).overlaps(1+h)
        s = 6*sum(one_position(i,shift) for i in range(4))
        costs.append(-to_bits*s)
    pair = sum(costs)
    assert 14400 < pair and pair < 14420
    multiplicity = [comb(256,k)*comb(256-k,k) for k in range(129)]
    assert all(b <= 2**16*a and a <= 2**13*b
               for a,b in zip(multiplicity,multiplicity[1:]))
    print("Volume-balanced pair suppression in bits:", pair)


if __name__ == "__main__":
    verify()
