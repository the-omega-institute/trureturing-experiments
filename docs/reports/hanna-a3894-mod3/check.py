#!/usr/bin/env python3
"""Exact finite checks of Hanna's A389472, A389474 and A389475 conjectures.
No network, third-party package, or floating-point arithmetic is used.
"""
from math import comb
from hashlib import sha256
import json

N = 403
SOURCE_474 = [1,1,4,14,104,488,3291,18141,127527,790427,5441663,35375651,244925743,1652908531,11518757017,79372305657,557324355555,3900077626815,27596212781442,195214162327620,1390657027410402,9919148868199578,71088535006091553,510378145935043827]
SOURCE_475 = [1,1,7,17,147,445,2199,9993,61837,281417,1515387,7844861,42920817,232527741,1297380639,7159441425,40375866078,226515725383,1286378694156,7313784351897,41879211301623,240239676368409,1384901325628119]
SOURCE_472 = [1,1,2,3,6,9,14,24,42,72,114,179,284,453,736,1213,2016,3362,5606,9318,15378,25121,40614,65216,104422,167118,267674,429297,689478,1109149,1788032,2889897,4683664,7610640,12394998,20225081,33050490,54068421,88519060,144984819,237508938]


def mul(a, b, degree, modulus=None):
    c = [0] * (degree + 1)
    for i, ai in enumerate(a[:degree+1]):
        if ai:
            for j, bj in enumerate(b[:degree+1-i]):
                if bj:
                    c[i+j] += ai * bj
        if modulus:
            c = [x % modulus for x in c]
    return c


def composition(a, h, degree):
    """Independent Horner evaluation, retaining all coefficients up to degree."""
    assert h[0] == 0
    order = next(i for i, v in enumerate(h) if v)
    c = [0] * (degree + 1)
    for k in range(degree // order, -1, -1):
        c = mul(c, h, degree)
        c[0] += a[k]
    return c


def a474(nmax):
    # Incremental coefficients of H^k; H has order two.
    a = [0] * (nmax + 2)
    a[1] = 1
    hp = [[0] * (nmax + 2) for _ in range((nmax + 1)//2 + 1)]
    h = hp[1]
    h[2] = 3
    for n in range(2, nmax+1):
        d = n+1
        for k in range(2, d//2+1):
            hp[k][d] = sum(h[j] * hp[k-1][d-j]
                           for j in range(2, d-2*(k-1)+1))
        nonlinear = sum(a[k] * hp[k][d] for k in range(2, d//2+1))
        square_known = sum(a[i] * a[d-i] for i in range(2, n))
        numerator = (2 if d == 3 else 0) - square_known + nonlinear
        assert numerator % 2 == 0, ('474-integrality', n)
        a[n] = numerator // 2
        h[d] = (2 if d == 3 else 0) + 2*(2*a[n]+square_known)
    return a[:nmax+1]


def a475(nmax):
    a = [0] * (nmax+2)
    b = [0] * (nmax+2)
    a[1] = b[2] = 1
    for n in range(2, nmax+1):
        d = n+1
        square_known = sum(a[i]*a[d-i] for i in range(2,n))
        cube = sum(a[i]*b[d-i] for i in range(1,d-1))
        inner = 3**(d//2)*a[d//2] if d % 2 == 0 else 0
        numerator = inner - square_known + 2*cube
        assert numerator % 2 == 0, ('475-integrality', n)
        a[n] = numerator // 2
        b[d] = 2*a[n]+square_known
    return a[:nmax+1]


def a472(nmax):
    a = [0,1,1]
    for m in range(3,nmax+1):
        d = m+2
        a.append(sum(a[k]*comb(k,d-2*k)
                     for k in range((d+2)//3,d//2+1)))
    return a


def frobenius(a):
    cube = mul(mul(a,a,N),a,N)
    sub = [a[d//3] if d%3 == 0 else 0 for d in range(N+1)]
    assert all((cube[d]-sub[d])%3 == 0 for d in range(N+1))
    return cube,sub


def main():
    a,b,c = a474(N),a475(N),a472(3*(N+1)-1)
    for name,seq,sample in [('A389474',a,SOURCE_474),('A389475',b,SOURCE_475),('A389472',c,SOURCE_472)]:
        assert seq[1:len(sample)+1] == sample, ('source-fidelity',name)
        digest = sha256(json.dumps(seq,separators=(',',':')).encode()).hexdigest()
        print(name,'source displayed terms matched =',len(sample),'; generated through index',len(seq)-1,'; coefficients sha256 =',digest)
    aa = mul(a,a,N)
    h = [2*x for x in aa]
    h[2] += 1
    h[3] += 2
    assert composition(a,h,N) == [3*x for x in aa]
    assert all(aa[d]%9 == 0 for d in range(4,N+1))
    print(f'A389474 (3): equation and square congruence at degrees 0..{N}; {N-3} requested coefficients at degrees 4..{N}')
    cube,sub = frobenius(a)
    p=q=[1]+[0]*N
    for k in range(1,N//3+2):
        p=mul(p,[x%3 for x in cube],N,3)
        q=mul(q,[x%3 for x in sub],N,3)
        assert p==q,('474 (4)',k)
    print(f'A389474 (4): degrees 0..{N} for k=1..{N//3+1}; k>={N//3+1} vanish in this truncation by order')
    bb=mul(b,b,N)
    bbb,subb=frobenius(b)
    lhs=[3**(d//2)*b[d//2] if d%2==0 else 0 for d in range(N+1)]
    rhs=[bb[d]-2*bbb[d]+(2 if d==2 else 0) for d in range(N+1)]
    assert lhs==rhs
    assert all((bb[d]-(1 if d==2 else 0)-2*subb[d])%3==0 for d in range(N+1))
    print(f'A389475 (4.a), (4.b): defining equation and both congruences checked at degrees 0..{N}')
    assert all(c[3*n-1]%3==0 for n in range(2,N+2))
    assert all(c[2*n-1]%2==0 for n in range(2,(len(c)+1)//2))
    # Direct binomial composition, checked independently of the recursion loop.
    for d in range(len(c)):
        comp=sum(c[k]*comb(k,d-2*k) for k in range((d+2)//3,d//2+1))
        assert comp==(c[d-2]+(1 if d==2 else 0) if d>=2 else 0)
    print(f'A389472 mod 3: {N} target indices n=2..{N+1}, coefficients through {3*(N+1)-1}; defining equation through degree {3*(N+1)-1}; existing parity line also checked')
    print('A389474 square coefficients, degrees 2..10:',','.join(map(str,aa[2:11])))
    print('ALL_CHECKS_OK')

if __name__=='__main__':
    main()
