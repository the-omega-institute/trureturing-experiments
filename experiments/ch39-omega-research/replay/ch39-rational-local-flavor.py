"""Finite Gaussian-rational shell and flavor readout, given two rational clock codes."""
from itertools import combinations, permutations, product
from math import isqrt
from pathlib import Path

from flint import fmpq as Q, fmpq_mat as Matrix

scale = 1 << 40
radius2 = 867
p = (1,5,29)
pole = (17,17,17)
portals = ((13,23,13),(23,13,13))
contacts = ((5,1,29),(11,11,25))
winner_codes = ((-1,-496481263766),(-1,489993214324))


def dot(v,w):
    return sum(x*y for x,y in zip(v,w))


def cross(v,w):
    return tuple(v[(i+1)%3]*w[(i+2)%3]-v[(i+2)%3]*w[(i+1)%3] for i in range(3))


def nearest(a,b):
    # Symmetric rounding ensures that matrix transposition commutes with rounding.
    return (1 if a >= 0 else -1)*((2*abs(a)+b)//(2*b))


def mat(rows,cols,entry):
    return Matrix([[entry(i,j) for j in range(cols)] for i in range(rows)])


def realify(a,b):
    r,c = a.nrows(),a.ncols()
    def entry(i,j):
        if (i<r)==(j<c):
            return a[i%r,j%c]
        return (-1 if i<r else 1)*b[i%r,j%c]
    return mat(2*r,2*c,entry)


def rounded(v):
    return Q(nearest(int(v.p)*scale,int(v.q)),scale)


def clock_pair(sign,n):
    den = scale*scale+n*n
    return (nearest(sign*scale*(scale*scale-n*n),den),
            nearest(sign*2*n*scale*scale,den))


def link(x,y,d,phases):
    a = scale*radius2*(radius2+dot(pole,x)+dot(x,y)+dot(y,pole))
    b = -isqrt(radius2*scale*scale)*dot(pole,cross(x,y))
    length = isqrt(a*a+b*b)
    u,v = nearest(scale*a,length),nearest(scale*b,length)
    if d == 200:
        u,v = -u,-v
    h,k = phases.get(d,(scale,0))
    return Q(-nearest(u*h-v*k,scale),scale), Q(-nearest(v*h+u*k,scale),scale)


def shell_response(codes):
    outer = sorted(set(permutations(p)))
    inner = sorted(set(permutations((13,13,23))))
    sites = outer+inner+sorted(v for v in product(range(30),repeat=3)
                               if dot(v,v)==radius2 and v not in outer+inner)
    size = len(sites)
    assert size == 19 and sites[0] == p
    re,im = Matrix(size,size),Matrix(size,size)
    phases = dict(zip((724,804),(clock_pair(*code) for code in codes)))
    for i,x in enumerate(sites):
        s = sum(x)
        re[i,i] = 20000*(s-35)*(s-49)**2 - 6 - Q(2*(s-35),7)
    for i,j in combinations(range(size),2):
        d = 2*(radius2-dot(sites[i],sites[j]))
        if d <= 804:
            re[i,i] += 1
            re[j,j] += 1
            r,t = link(sites[i],sites[j],d,phases)
            re[i,j],im[i,j],re[j,i],im[j,i] = r,t,r,-t
    internal = realify(mat(18,18,lambda i,j:re[i+1,j+1]), mat(18,18,lambda i,j:im[i+1,j+1]))
    a,b = Matrix(3,18),Matrix(3,18)
    for end in portals:
        j = sites.index(end)
        n = cross(p,end)
        length = isqrt(scale*scale*dot(n,n))
        for i in range(3):
            unit = Q(n[i]*scale,length)
            a[i,j-1], b[i,j-1] = unit*re[0,j],unit*im[0,j]
    f = realify(a,b)
    native = -f*internal.solve(f.transpose())
    return mat(6,6,lambda i,j:rounded(native[i,j]))


def projector(n):
    return mat(3,3,lambda i,j:Q(n[i]*n[j],dot(n,n)))


def flavor_operation(codes):
    eye = mat(3,3,lambda i,j:int(i==j))
    zero = Matrix(3,3)
    plane = realify(eye-projector(p),zero)
    q,r = [realify(projector(cross(p,b)),zero) for b in contacts]
    g,h,k,alpha = [Q(n,2048) for n in (1486,5202,5069,6322)]
    action = g*plane*shell_response(codes)*plane+h*q+k*r+alpha*plane
    assert action == action.transpose()
    unit = realify(eye,zero)
    z = realify(zero,eye)*action/16
    even = 120*unit+12*z**2
    odd = 60*z+z**3
    w = (even-odd).solve(even+odd)**16
    assert w*w.transpose() == unit
    return w


def measurements(w):
    def prob(i,j):
        return w[i,j]**2+w[i+3,j]**2
    u2 = prob(0,2)
    c2 = 1-u2
    s2,t2 = prob(0,1)/c2,prob(1,2)/c2
    # J = Im((V00 conjugate(V10)) conjugate(V01 conjugate(V11))).
    def pair(j):
        a,b,c,d = w[0,j],w[3,j],w[1,j],w[4,j]
        return a*c+b*d, b*c-a*d
    ar,ai = pair(0)
    br,bi = pair(1)
    jcp = ai*br-ar*bi
    cos_term = (prob(1,0)-s2*(1-t2)-(1-s2)*t2*u2)/2
    return (s2,t2,u2), jcp, cos_term


def records():
    table = Path(__file__).with_name("ch39-flavor-clock-codes.tsv").read_text().splitlines()
    assert table[0].split() == ["index","sign724","code724","sign804","code804"]
    assert len(table) == 1279
    for idx,line in enumerate(table[1:],10815):
        n,s,a,t,b = map(int,line.split())
        assert n == idx and s in (-1,1) and t in (-1,1) and max(abs(a),abs(b)) <= scale
        yield idx,((s,a),(t,b))


def verify():
    half = Q(55,1000000)
    centers = (Q(22501,100000),Q(4183,100000),Q(3732,1000000))
    winners = []
    for idx,codes in records():
        squares,jcp,cos_term = measurements(flavor_operation(codes))
        if any(z < (c-half)**2 or z > (c+half)**2 for c,z in zip(centers,squares)):
            continue
        if jcp <= 0 or cos_term <= 0:
            continue
        tan_delta = jcp/((1-squares[2])*cos_term)
        if tan_delta < Q(2181,1000) or tan_delta > Q(2253,1000):
            continue
        assert Q(2182,1000) < tan_delta < Q(2252,1000)
        winners.append(idx)
        assert codes == winner_codes
        print("squared sines and J, rounded at 40 bits:",
              *(nearest(int(v.p)*scale,int(v.q)) for v in (*squares,jcp)))
    assert winners == [11431]
    print("all-rational 1278-code scan:", winners)


if __name__ == "__main__":
    verify()
