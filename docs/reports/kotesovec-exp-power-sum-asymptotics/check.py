#!/usr/bin/env python3
"""Exact checks of seven Hanna series and Kotesovec asymptotic squeezes.

OEIS term data: Paul D. Hanna / OEIS Foundation, OEIS End-User License.
Uses integers and fractions.Fraction only; finite checks do not prove limits.
"""
from fractions import Fraction as F
from math import comb, factorial

EXPECTED = {"A156170": [1, 1, 3, 10, 41, 219, 1602, 16635, 247171, 5242108, 157390565, 6663089873, 396778864166, 33200932308437, 3906922702271961], "A276751": [1, 1, 3, 14, 111, 1813, 57846, 3941129, 515554887, 139563384274, 73929755773659, 78682910542834037, 169524995438153307498, 712160156293232925362965, 6241130803695426404771763891], "A276752": [1, 1, 5, 30, 327, 7085, 307280, 28472653, 5000661017, 1886425568702, 1331753751874235, 2008313162512681569, 5765904212733638946976, 34525801618218187545094977, 406111805399407205212602871837], "A292500": [1, 1, 4, 18, 122, 1382, 26992, 967860, 59207134, 6539607238, 1225903048760, 407719392472476, 233686070341415140, 233030334505100451484, 407716349332865096406960], "A292502": [1, 7, 43, 399, 6091, 151255, 6550307, 465127199, 58293976795, 12191724780647, 4471204259257363, 2799295142330495151, 3026340345288168023883, 5704756586858875194533367, 18287793731664040419412785283], "A292501": [1, 1, 3, 13, 91, 1119, 23235, 879361, 55447631, 6274018595, 1192773105789, 400761393446831, 231147252957096671, 231434829013884972151, 406000810484101907916927], "A156100": [1, 1, 3, 7, 25, 113, 741, 7181, 101139, 2089283, 61683087, 2600572391, 156100460443, 13231060891179, 1594932996895155]}
LIMIT = 40


def mul(a, b):
    n = len(a)
    out = [0] * n
    for i, x in enumerate(a):
        if x:
            for j in range(n - i):
                if b[j]:
                    out[i + j] += x * b[j]
    return out


def power(a, m):
    out = [1] + [0] * (len(a) - 1)
    for _ in range(m):
        out = mul(out, a)
    return out


def inner(anum, m, n):
    if anum == 'A156100':
        return [0, 1, 2**m] + [0] * (n - 2)
    exponent = 2*m-1 if anum == 'A276751' else 2*m if anum == 'A276752' else m
    odd = anum in ('A292500', 'A292501', 'A292502')
    return [0] + [(2*k-1 if odd else k)**exponent for k in range(1, n+1)]


def exponential(b):
    a = [F(1)]
    for n in range(1, len(b)):
        a.append(sum((k*b[k]*a[n-k] for k in range(1, n+1)), F(0))/n)
    return a


def calculate(anum):
    n = LIMIT
    b = [F(0)] * (n+1)
    for m in range(1, n+1):
        original = inner(anum, m, n)
        p = power(original, m)
        if anum == 'A292501':
            # Direct signed source expression versus the positive first difference.
            one_minus_x = [1, -1] + [0]*(n-1)
            signed = mul(p, power(one_minus_x, m))
            delta = [0]+[original[k]-original[k-1] for k in range(1,n+1)]
            assert all(v >= 0 for v in delta)
            p = power(delta, m)
            assert p == signed
        if anum == 'A156100':
            assert p == [comb(m,k-m)*2**(m*(k-m)) if m <= k <= 2*m else 0
                         for k in range(n+1)]
        for k in range(n+1):
            b[k] += F(p[k],m)
    a = [k*b[k] for k in range(n+1)] if anum == 'A292502' else exponential(b)
    start = 1 if anum == 'A292502' else 0
    assert a[start:start+15] == EXPECTED[anum], (anum, 'OEIS mismatch')
    assert all(x.denominator == 1 for x in a)
    assert all(x > 0 for x in a[1:])
    if anum != 'A292502':
        # Independent finite exponential expansion, not the derivative recurrence.
        short = b[:16]
        direct = [F(0)] * 16
        for r in range(16):
            p = power(short,r)
            direct = [x + F(y,factorial(r)) for x,y in zip(direct,p)]
        assert a[:16] == direct
    return b,a


def lower(anum,n):
    m,t=(n+1)//2,n//2
    if anum in ('A156170','A156100'):return F(2**(m*t),m)
    if anum == 'A276751':return F(2**((2*m-1)*t),m)
    if anum == 'A276752':return F(2**(2*m*t),m)
    if anum == 'A292501':return F((3**m-1)**t,m)
    if anum == 'A292502':return F(n*3**(m*t),m)
    return F(3**(m*t),m)


def check_maxima(anum):
    # Independent max-product dynamic program over every composition, N <= 20.
    for m in range(1,21):
        weights = inner(anum,m,20)
        if anum == 'A292501':weights=[0]+[weights[k]-weights[k-1] for k in range(1,21)]
        best=[1]+[0]*20
        for _ in range(m):
            best=[0]+[max((best[n-k]*weights[k] for k in range(1,n+1)),default=0)
                      for n in range(1,21)]
        base=3 if anum in ('A292500','A292501','A292502') else 2
        scale=2 if anum in ('A276751','A276752') else 1
        for n in range(m,21):
            assert best[n] <= base**(scale*m*(n-m))
            assert m*(n-m) <= n*n//4


def main():
    for anum in EXPECTED:
        b,a=calculate(anum)
        base=3 if anum in ('A292500','A292501','A292502') else 2
        scale=2 if anum in ('A276751','A276752') else 1
        check_maxima(anum)
        for n in range(1,LIMIT+1):
            q=base**(scale*(n*n//4))
            assert 0 < b[n] <= 2**(n-1)*q
            upper = n*2**(n-1)*q if anum == 'A292502' else 4**n*q
            assert lower(anum,n) <= a[n] <= upper, (anum,n)
            if anum == 'A292501':
                m,t=(n+1)//2,n//2
                assert lower(anum,n) >= F(2,3)**t*F(3**(m*t),m)
        print(anum, 'FIDELITY first15 =', EXPECTED[anum])
        print('  Exact lower <= a(N) <= upper and B bounds: every 1 <= N <=',LIMIT)
        print('  Max-product bounds: every 1 <= m <= N <= 20')
        for n in (1,2,3,5,10,20,30,40):
            print('  N =',n,'lower =',lower(anum,n),'PASS')
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
