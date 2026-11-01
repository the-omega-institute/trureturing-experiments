#!/usr/bin/env python3
"""Exact checks for Hanna's A352409 COMMENTS (4); no assertion of a uniform proof.

Source: https://oeis.org/search?q=id:A352409&fmt=text (revision 17).
Companion: https://oeis.org/search?q=id:A352408&fmt=text (revision 12).
Published finite prefix attributed to Paul D. Hanna; OEIS End-User License:
https://oeis.org/LICENSE . All computation uses Python integer/Fraction arithmetic.
Run: python3 check.py. Standard library only; no network access is needed.
"""
from fractions import Fraction
from math import comb, gcd
from random import Random

SOURCE_PREFIX = (1, 3, 45, 1267, 51597, 2761539, 182885885, 14415019395, 1316237331069, 136512958750979, 15842506286290173, 2033176597680449283, 285833727841312233725, 43677225803116362273795, 7207197437612731825348605, 1277141936892060488486787075)
SOURCE_ROWS = (
 (1, 0, 0, 0, 0, 0, 0),
 (1, 3, 0, -14, -153, -4059, -162214),
 (1, 10, 35, 0, -825, -17758, -642015),
 (1, 21, 189, 847, 0, -55818, -1835218),
 (1, 36, 594, 5772, 32715, 0, -4524660),
 (1, 55, 1430, 23100, 252450, 1762706, 0),
 (1, 78, 2925, 69836, 1179672, 14597856, 122423756),
)
COMPANION_PREFIX = (0, 3, 105, 4521, 247509, 16743276)


def mul(f, g, d):
    out = [0] * (d + 1)
    for i, a in enumerate(f[:d+1]):
        for j, b in enumerate(g[:d+1-i]):
            out[i+j] += a*b
    return out


def power(f, e, d):
    out = [1] + [0]*d
    while e:
        if e % 2:
            out = mul(out, f, d)
        e //= 2
        if e:
            f = mul(f, f, d)
    return out


def power_by_derivative(f, e, d):
    # Independent coefficient recurrence: f*(f^e)' = e*f'*f^e.
    assert f[0] == 1
    out = [1]
    for j in range(1, d+1):
        numerator = sum(((e+1)*i-j)*f[i]*out[j-i]
                        for i in range(1, min(j, len(f)-1)+1))
        q, r = divmod(numerator, j)
        assert r == 0
        out.append(q)
    return out


def normalized(f, M, m, independent=False):
    assert m > 0 and M > 0 and gcd(M,m) == 1 and f[0] == 1
    calc = power_by_derivative if independent else power
    q, r = divmod(calc(f, M, m)[m], M)
    assert r == 0, (M,m,r)
    return q


def source_series(count):
    a = []
    for n in range(1, count+1):
        m, M = n+1, n*(2*n+1)
        h = [1,1] + [-v for v in a] + [0]
        value = normalized(h,M,m)
        assert value == normalized(h,M,m,independent=True)
        a.append(value)
        h[m] = -value
        assert power(h,M,m)[m] == 0
    return a


def main():
    a = source_series(40)
    assert tuple(a[:len(SOURCE_PREFIX)]) == SOURCE_PREFIX
    for n, row in enumerate(SOURCE_ROWS):
        assert tuple(power([1,1]+[-v for v in a],n*(2*n+1),6)) == row
    assert all(v % 2 == 1 for v in a)
    print('FIDELITY: reproduced all 16 displayed A352409 terms and 7 displayed rows (degrees 0..6).')
    print('SOURCE: independently generated a_0..a_39; both algorithms agree; all 40 are odd; all 40 defining diagonals vanish.')
    print('FIRST_8 =', a[:8])

    residuals=[]
    companion=[]
    for n in range(1,101):
        m,M=n+1,n*(2*n+1)
        assert M % m == 1 and gcd(M,m) == 1
        h=[1,1]+[-1]*(m-1)
        direct=normalized(h,M,m)
        assert direct == normalized(h,M,m,independent=True)
        terms=[Fraction((-2)**k*comb(M,k)*comb(M+n-2*k,m-2*k),M)
               for k in range(m//2+1)]
        assert sum(terms) == direct
        assert terms[0] == Fraction(2*comb(2*n*(n+1)-1,n),2*n+1)
        assert all(t.numerator % 2 == 0 and t.denominator % 2 == 1 for t in terms)
        assert direct % 2 == 0
        residuals.append(direct)
        h[m] += 1
        assert normalized(h,M,m) == direct+1
        c3=[1,1]+[-comb(2*(j//3),j//3)//(j//3+1) if j%3==0 else 0
                     for j in range(m-1)]
        r3=normalized(c3,M,m)
        assert r3 == normalized(c3,M,m,independent=True)
        assert r3%3 == 0
        companion.append(r3)
    assert tuple(companion[:6]) == COMPANION_PREFIX
    print('PARITY_REFERENCE: exact expansion, leading binomial identity, termwise even localization, and residual evenness passed for n=1..100.')
    print('R_1..R_6 =', residuals[:6])
    print('DEGREE_m_UPDATE: T(H+x^m)=T(H)+1 passed for n=1..100.')
    print('MOD3_TRIAGE: reproduced A352408 terms 1..6 and tested residual divisibility for n=1..100; finite evidence only.')
    for j,value in enumerate(a):
        target=comb(2*(j//3),j//3)//(j//3+1) if j%3==0 else 0
        assert (value-target)%3 == 0
    print('MOD3_TRIAGE: all 40 generated A352409 coefficients satisfy COMMENTS (1)..(3); no uniform mod3 theorem asserted.')

    rng=Random(352409)
    trials=0
    for p in (2,3):
        for m in range(1,9):
            for M in range(1,19):
                if gcd(M,m)!=1:
                    continue
                g=[1]+[rng.randrange(-3,4) for _ in range(m)]
                f=[1]+[g[i]+p*rng.randrange(-3,4) for i in range(1,m+1)]
                assert (normalized(f,M,m)-normalized(g,M,m))%p == 0
                trials+=1
    for p in (2,3):
        for k in range(1,257):
            q=Fraction(p**k,k)
            assert q.numerator%p==0 and q.denominator%p!=0
    print(f'NORMALIZED_TRANSFER: {trials} exact random pairs (seed 352409) passed modulo 2 and 3; p^k/k localization checked for k=1..256.')
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
