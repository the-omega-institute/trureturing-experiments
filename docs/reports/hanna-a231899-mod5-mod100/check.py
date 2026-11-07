#!/usr/bin/env python3
"""Exact finite checks for Hanna's OEIS A231899 modulo-5/modulo-100 conjectures.

Python standard library only. This does not prove an unbounded assertion.
Source: https://oeis.org/A231899, approved revision 19, and its b-file.
OEIS data remain under https://oeis.org/LICENSE; new code: MIT.
"""
from fractions import Fraction as Q
from math import comb, factorial

N = 40
DATA = (1,1,1,7,106,2590,89500,4061830,230971450,15931539370,
        1301850919000,123784442479000,13503177242345500,1670758163868804100,
        232275860883184444000,35996619826480459259200,6176482540368771725898700,
        1166516670494844414929980300,241253958437270751756621102400)


def egf_mul(a, b, n):
    return [sum(comb(k,j)*a[j]*b[k-j]
                for j in range(max(0,k-len(b)+1),min(k+1,len(a))))
            for k in range(n+1)]


def egf_pow(a, e, n):
    out = [1]+[0]*n
    while e:
        if e % 2:
            out = egf_mul(out,a,n)
        e //= 2
        if e:
            a = egf_mul(a,a,n)
    return out


def bell_sequence(nmax):
    """Triangular update from A'=A(x(A')^3), using partial Bell polynomials."""
    a = [1]
    for n in range(nmax):
        cube = egf_pow(a[1:],3,max(0,n-1))
        h = [0]+[j*cube[j-1] for j in range(1,n+1)]
        bell = [[0]*(n+1) for _ in range(n+1)]
        bell[0][0] = 1
        for t in range(1,n+1):
            for k in range(1,t+1):
                bell[t][k] = sum(comb(t-1,j-1)*h[j]*bell[t-j][k-1]
                                 for j in range(1,t-k+2))
        a.append(sum(a[k]*bell[n][k] for k in range(n+1)))
    return a


def lagrange_sequence(nmax):
    """Independent source formula: a_n=E_(n-1)(A^(3n-2))/(3n-2)."""
    a = [1]
    for n in range(1,nmax+1):
        v = egf_pow(a,3*n-2,n-1)[n-1]
        q,r = divmod(v,3*n-2)
        assert r == 0, (n,r)
        a.append(q)
    return a


def mul(a,b,limit=None):
    size = len(a)+len(b)-1
    if limit is not None:
        size = min(size,limit+1)
    out = [Q(0)]*size
    for i,x in enumerate(a):
        for j,y in enumerate(b[:max(0,size-i)]):
            out[i+j] += x*y
    return out


def compose(f,h,limit=None):
    out = [Q(0)]
    for c in reversed(f):
        out = mul(out,h,limit)
        out[0] += c
    while len(out)>1 and out[-1]==0:
        out.pop()
    return out


def derivative(f):
    return [j*f[j] for j in range(1,len(f))]


def residual(f):
    b = derivative(f)
    h = [Q(0)]+mul(mul(b,b),b)
    comp = compose(f,h)
    return [(b[n] if n<len(b) else Q(0))-comp[n] for n in range(len(comp))]


def vp(q,p):
    assert q != 0
    num,den = abs(q.numerator),q.denominator
    v = 0
    while num % p == 0:
        num //= p
        v += 1
    while den % p == 0:
        den //= p
        v -= 1
    return v


def main():
    a = bell_sequence(N)
    assert a == lagrange_sequence(N)
    assert tuple(a[:len(DATA)]) == DATA
    assert tuple(a) == BFILE_0_40
    assert all(x>=0 for x in a)
    assert all(a[n] % 5 == 0 for n in range(5,N+1))
    assert all(a[n] % 100 == 0 for n in range(10,N+1))
    print('Fidelity: all 19 DATA terms and b-file indices 0..40 agree.')
    print('Independent Bell and source Lagrange recurrences agree at 0..40.')
    print('Congruences: mod 5 at 5..40; mod 100 at 10..40.')
    print('Initial a_0..a_10:',a[:11])
    expected = [-1301850919000,-83427063990000,-7576062984979500,
                -797207176461504600,-92429727810268992000,
                -11568255207014755048800,-1546549868706068453584500]
    for degree,modulus,res_degree in [(4,5,40),(9,100,225)]:
        f = [Q(a[j],factorial(j)) for j in range(degree+1)]
        r = residual(f)
        assert len(r)-1 == res_degree
        normalized = [q*factorial(n) for n,q in enumerate(r)]
        assert all(q.denominator==1 and q.numerator%modulus==0 for q in normalized)
        assert all(q.denominator%5!=0 for q in f+derivative(f)+r)
        print(f'F_{degree}: coefficients =',','.join(map(str,f)))
        print(f'R_{degree}: degree {res_degree}; all {len(r)} normalized coefficients integral and divisible by {modulus}.')
        if degree == 4:
            assert normalized[:5] == [0,0,0,0,-2590]
        else:
            assert normalized[:9] == [0]*9
            assert normalized[9:16] == expected
            print('R_9 normalized coefficients at 9..15:',','.join(str(q.numerator) for q in normalized[9:16]))
            vals_f = [vp(q,2) for q in f]
            vals_b = [vp(q,2) for q in derivative(f)]
            assert vals_f == [0,0,-1,-1,-2,-2,-2,-3,-6,-6]
            assert vals_b == [0,0,-1,0,-2,-1,-3,-3,-6]
            assert all(v>=-(3*j//4) for j,v in enumerate(vals_f))
            assert all(v>=-(3*j//4) for j,v in enumerate(vals_b))
            assert all(q==0 or vp(q,2)>=-(3*j//4) for j,q in enumerate(r))
            print('v2(F_9 coefficients):',vals_f)
            print('v2(F_9 derivative coefficients):',vals_b)
    for n in range(16,10001):
        bound = n//2+n//4+n//8+n//16-3*n//4
        assert bound>=2
    print('Factorial lower-bound inequality checked for n=16..10000 (finite supplement).')
    # Check the supplied exclusions of two general congruence theorems.
    f = [Q(a[j],factorial(j)) for j in range(18)]
    log1p = [Q(0)]+[Q((-1)**(j+1),j) for j in range(1,5)]
    assert compose(f[:5],log1p,4)[4] == Q(23,8)
    # x A'/A = sum_m w_m x^m, computed without floating point or logarithms.
    inv = [Q(1)]
    for n in range(1,18):
        inv.append(-sum(f[j]*inv[n-j] for j in range(1,n+1)))
    w = [Q(0)]+mul(derivative(f),inv,16)
    assert w[4] == Q(27,2)
    assert w[17] == Q(114101705219927291073,2242240)
    assert w[17].denominator%5 == 0
    print('General-theorem exclusions: [y^4]A(log(1+y))=23/8; w_4=27/2; w_17='+str(w[17])+'.')
    print('These finite checks do not establish the unbounded conjectures.')
    print('ALL_CHECKS_OK')


# First 41 entries of https://oeis.org/A231899/b231899.txt, retained for offline reproduction.
BFILE_0_40 = (1, 1, 1, 7, 106, 2590, 89500, 4061830, 230971450, 15931539370, 1301850919000, 123784442479000, 13503177242345500, 1670758163868804100, 232275860883184444000, 35996619826480459259200, 6176482540368771725898700, 1166516670494844414929980300, 241253958437270751756621102400, 54389029505329397128287295934200, 13311650176211618048628173480676400, 3524050074283492777996348123378158400, 1005762460795500781061289014834083951600, 308511375293925305238589078083836048018800, 101427691966601629544138851621666113168621400, 35647772293993918983858541638057410226343355800, 13361700359547512381940828655749854256960508580800, 5329422114667013306164035065267679838790745750713200, 2257292500714313905873293382131843357349921038940471600, 1013312078570952008585465821942403605036055278500091478800, 481234799209362473390190433721911293085892568517417159021600, 241372393829045254103814517105449042486135035030996568622489600, 127654651970801311301783622622248733613614034726572852039026612400, 71079692943490996569941128032489134768807987065710649852682698433200, 41609458040145283290896135380672334438817226897037529540585736270441600, 25573340102342550866069982754458935063619434599377317524631830838916863200, 16480632682528642247925258642524378508124762172803325743968778766423135507200, 11123002770014142804444016133606542930273039365033753140560609417332099603476800, 7852882772519565357118628678869969207627446193828482085669483653559479787303614400, 5793177059071471746827479139494036032666008795699168475213593261805752071598943707200, 4460973547909610866201517278047075815687204208589965908541819843437086865421188945970400)

if __name__ == '__main__':
    main()
