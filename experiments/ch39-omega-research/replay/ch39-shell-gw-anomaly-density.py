"""Haar area of the anomaly-compatible two-clock set using certified Bernstein signs."""
from collections import deque
from math import comb, floor, log2
from pathlib import Path
from runpy import run_path
from flint import arb, arb_mat, acb_mat, ctx
import numpy as np


model=run_path(str(Path(__file__).with_name("ch39-shell-gw-rh-flavor.py")))
ctx.prec=200
degree=4
width=1024


def coefficients(base):
    points=[arb(i)/degree for i in range(degree+1)]
    angles=[(2*x-1).acos() for x in points]
    collocation=arb_mat([[comb(degree,j)*x**j*(1-x)**(degree-j) for j in range(degree+1)] for x in points])
    vals=[]
    for p in angles:
        row=[]
        for q in angles:
            det=acb_mat(model["family_matrix"](base,p,q)).det()
            assert det.imag.contains(0)
            row.append(det.real)
        vals.append(row)
    inv=collocation.inv()
    poly=inv*arb_mat(vals)*inv.transpose()
    return [[poly[i,j] for j in range(5)] for i in range(5)]


def halve(row):
    left,right=[],[]
    while len(row)>1:
        left.append(row[0]);right.append(row[-1])
        row=[(a+b)/2 for a,b in zip(row,row[1:])]
    return left+row,row+right[::-1]


def split(poly):
    columns=list(zip(*poly))
    xhalf=[halve(col) for col in columns]
    pieces=[]
    for i in (0,1):
        halfrows=list(zip(*(c[i] for c in xhalf)))
        yhalf=[halve(row) for row in halfrows]
        pieces.extend([[row[j] for row in yhalf] for j in (0,1)])
    return pieces


def regular(poly):
    return all(v>0 for row in poly for v in row) or all(v<0 for row in poly for v in row)


def estimate(blocks,polynomials):
    cells=bytearray(width**2)
    def visit(x,y,side,polys):
        undecided=[a for a in polys if not regular(a)]
        if not undecided:
            for i in range(x,x+side):
                cells[i*width+y:i*width+y+side]=bytes([1])*side
        elif side>1:
            h=side//2
            for k,children in enumerate(zip(*(split(a) for a in undecided))):
                visit(x+(k//2)*h,y+(k%2)*h,h,children)
    visit(0,0,width,polynomials)
    def theta(i):
        return (arb(2*i)/width-1).acos()
    weights=[(theta(i)-theta(i+1))/arb.pi() for i in range(width)]
    yes=arb(0);no=arb(0);unknown=arb(0);components=0
    for i in range(width):
        for j in range(width):
            start=i*width+j
            if cells[start]==0:
                unknown+=weights[i]*weights[j]
                continue
            if cells[start]!=1:
                continue
            components+=1
            queue=deque([start]); cells[start]=2
            chosen=start
            mass=arb(0)
            while queue:
                cell=queue.popleft()
                x,y=divmod(cell,width)
                mass+=weights[x]*weights[y]
                for a,b in ((x-1,y),(x+1,y),(x,y-1),(x,y+1)):
                    nb=a*width+b
                    if 0<=a<width and 0<=b<width and cells[nb]==1:
                        cells[nb]=2
                        queue.append(nb)
            x,y=divmod(chosen,width)
            p=(arb(2*x+1)/width-1).acos()
            q=(arb(2*y+1)/width-1).acos()
            kernels=[model["family_matrix"](b,p,q) for b in blocks]
            indices=tuple(18-model["model"]["count"](k,arb(0)) for k in kernels)
            if indices==(3,-3,0)[:len(blocks)]:
                for k,index in zip(kernels,indices):
                    guess=np.linalg.eigvalsh([[complex(x) for x in row] for row in k])
                    gap=arb(2)**(floor(log2(min(abs(guess))))-1)
                    assert model["model"]["count"](k,-gap)==model["model"]["count"](k,gap)
                    model["opposite_modes"](k,index,gap)
                yes+=mass
            else:
                no+=mass
    assert (yes+no+unknown).contains(1)
    return yes,1-no,components


if __name__=="__main__":
    blocks=[model["model"]["make"](sign,True) for sign in (1,-1,0)]
    polynomials=[coefficients(b) for b in blocks]
    for count in (2,3):
        lo,hi,components=estimate(blocks[:count],polynomials[:count])
        print(count,"tested sectors; certified density",lo.lower().str(12),hi.upper().str(12),
              "connected certified regions",components)
