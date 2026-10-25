#!/usr/bin/env python3
"""Exact checks of five Peter Bala OEIS phi(k)^2 conjectures.

Run with Python 3, standard library only. All evidence uses integers or Fraction.
The finite checks support a preregistration; they do not prove eventual periodicity.
"""
from fractions import Fraction as F
from math import factorial

PREFIX = {'A281180': [1, 4, 88, 4672, 454144, 70084096, 15728822272, 4836914249728, 1952137912385536, 1000749157519458304, 635146072839001735168, 488855521088102855606272, 448599416591747486039670784, 483861305506660094099058589696, 606050665000453965359938841608192, 872366179652871528356910686198038528, 1430068361869553198039835379199635357696, 2648687881942689612933392158083076801429504, 5503854158077547090902251582359116752300802048], 'A281181': [1, 1, 13, 493, 37369, 4732249, 901188997, 240798388357, 85948640603761, 39504564917358001, 22726779729476308093, 15998009117983994065693, 13526765851190230940840809, 13528070218935445806530640649, 15795819619923464298050697616117, 21294937666865806704402646632389557, 32828500597549179599563478551377297121, 57385924456400269824204023290894357442401, 112904615348383588847189789579363784912180973], 'A281183': [1, 2, 32, 1376, 114176, 15519488, 3132551168, 879422726144, 327670676455424, 156439068819587072, 93116847688811282432, 67602541384815095054336, 58796336342280763841970176, 60351125684887424790500999168, 72187248798124538021926003539968, 99529442030183464236437157900713984, 156697512616609083360755035696397287424], 'A281184': [1, 3, 57, 2739, 246801, 35822307, 7636142793, 2246286827091, 871869519033249, 431649452286233283, 265466419357802436057, 198541440131880248161779, 177448471205103040365902001, 186781461066456338787698757027, 228695537454759099917373077023593, 322272887805877963568678968978370451, 517868815187736150011294497645677002049], 'A281440': [1, 1, 29, 2605, 478745, 148838425, 70130095925, 46612385237125, 41546401457128625, 47826888286436568625, 69072143074283849778125, 122288686408468623492188125, 260460302659740930298833415625, 657013212464002677825677944215625, 1937092871632377472727255256840753125, 6600357306119497177404312427298619203125, 25738551995192677896309032835665731654390625]}
MEMBERS = {"A281181": (3, 1), "A281183": (3, 2),
           "A281184": (3, 3), "A281440": (7, 1)}


def mul(a, b, degree):
    out = [F(0)] * (degree + 1)
    for i, x in enumerate(a[:degree + 1]):
        for j, y in enumerate(b[:degree + 1 - i]):
            out[i + j] += x * y
    return out


def power(a, e, degree):
    out = [F(1)]
    while e:
        if e & 1:
            out = mul(out, a, degree)
        a = mul(a, a, degree)
        e //= 2
    return out


def source_series(q, degree):
    # Literal coefficient recursion for C'=C^q S, S'=C^(q+1), C(0)=1,S(0)=0.
    c, s = [F(1)], [F(0)]
    for n in range(degree):
        cq = power(c, q, n)
        cp = mul(cq, s, n)[n] / (n + 1)
        sp = mul(cq, c, n)[n] / (n + 1)
        c.append(cp)
        s.append(sp)
    return c, s


def d_integer(poly):
    # D=(1+y^2)^2 d/dy on Z[y].
    out = [0] * (len(poly) + 3)
    for j, a in enumerate(poly):
        if j:
            for shift, factor in [(0, 1), (2, 2), (4, 1)]:
                out[j - 1 + shift] += j * a * factor
    return out


def e_integer(poly, q):
    # E=D^2 on Z[c], using s^2=c^2-1; D c=c^q s, D s=c^(q+1).
    out = [0] * (len(poly) + 2 * q)
    for j, a in enumerate(poly):
        out[j + 2 * q] += a * j * (j + q)
        out[j + 2 * q - 2] -= a * j * (j + q - 1)
    return out


def fidelity():
    # Compare two independent constructions, then all displayed OEIS terms.
    c3, s3 = source_series(3, 38)
    c7, _ = source_series(7, 36)
    d, odd = [0, 1], []
    for m in range(1, 2 * len(PREFIX['A281180'])):
        d = d_integer(d)
        if m % 2:
            odd.append(d[0])
    assert odd == PREFIX['A281180']
    for n, a in enumerate(odd, 1):
        assert factorial(2 * n - 1) * s3[2 * n - 1] == a
    for name, (q, exponent) in MEMBERS.items():
        p = [0] * exponent + [1]
        vals = []
        series = power(c3 if q == 3 else c7, exponent, 36)
        for n, want in enumerate(PREFIX[name]):
            got = sum(p)
            assert got == want, (name, n, got, want)
            assert factorial(2 * n) * series[2 * n] == got
            vals.append(got)
            p = e_integer(p, q)
        print('FIDELITY', name, len(vals), 'displayed terms; rational EGF cross-check OK')
    # Test the exact hypotheses of the screened integer-logarithmic-weight theorems.
    for name, source, exponent, index, expected in [
        ('A281180 derivative Sprime', c3, 4, 4, F(20,3)),
        ('A281181', c3, 1, 4, F(5,3)),
        ('A281183', c3, 2, 4, F(10,3)),
        ('A281184', c3, 3, 6, F(41,5)),
        ('A281440', c7, 1, 4, F(13,3))]:
        f = power(source, exponent, 6)
        w4 = 4 * (f[4] - f[2]**2 / 2)
        w6 = 6 * (f[6] - f[2]*f[4] + f[2]**3 / 3)
        weight = w4 if index == 4 else w6
        assert weight == expected and weight.denominator != 1
        print('LOG_WEIGHT', name, 'w'+str(index)+'='+str(weight), 'nonintegral')
    block = [1,4,4,3,5,5,2,1,1,6,3,3,4,2,2,5,6,6]
    assert [a % 7 for a in odd[:18]] == block
    assert all(any(block[i] != block[(i+t) % 18] for i in range(18))
               for t in range(1, 18))
    print('FIDELITY A281180', len(odd), 'displayed terms; rational EGF cross-check OK')
    print('FIDELITY A281180 modulo 7: displayed 18-term block; least cyclic period 18')
    assert [a % 13 for a in PREFIX['A281440'][1:13]] == [1,3,5,7,8,6,6,5,4,3,9,10]
    print('SOURCE_INDEX A281440 modulo 13: source example begins at n=1; offset-0 prefix is [1,1,3]')


def factors(n):
    out = []
    p = 2
    while p * p <= n:
        if n % p == 0:
            r = 0
            while n % p == 0:
                n //= p
                r += 1
            out.append((p, r))
        p += 1
    if n > 1:
        out.append((n, 1))
    return out


def phi(n):
    out = n
    for p, _ in factors(n):
        out = out // p * (p - 1)
    return out


def modular_s(modulus, end):
    # The ideal (y^k) is D-stable modulo k: D(y^k)=0 there.
    if modulus == 1:
        return [0] * (end + 1)
    p = [0] * modulus
    p[1] = 1
    out = [0] * (end + 1)
    for m in range(1, 2 * end):
        t = [0] * modulus
        for j, a in enumerate(p):
            if a and j:
                for shift, factor in [(0, 1), (2, 2), (4, 1)]:
                    dest = j - 1 + shift
                    if dest < modulus:
                        t[dest] += j * a * factor
        p = [a % modulus for a in t]
        if m % 2:
            out[(m + 1) // 2] = p[0]
    return out


def modular_c(modulus, end, q, exponent):
    # The ideal (c^k-1) is E-stable modulo k, so evaluation at c=1 is preserved.
    if modulus == 1:
        return [0] * (end + 1)
    p = [0] * modulus
    p[exponent % modulus] = 1
    out = []
    for _ in range(end + 1):
        out.append(sum(p) % modulus)
        t = [0] * modulus
        for j, a in enumerate(p):
            t[(j + 2 * q) % modulus] += a * j * (j + q)
            t[(j + 2 * q - 2) % modulus] -= a * j * (j + q - 1)
        p = [a % modulus for a in t]
    return out


def base_operator_check():
    # Exact characteristic-p identities on F_p[y]/(y^p), for both polynomial flows.
    for p in [2,3,5,7,11,13,17,19]:
        for power_exponent in [2,4]:
            from math import comb
            pp = [(2*i, comb(power_exponent, i)) for i in range(power_exponent+1)]
            def apply(v):
                w = [0] * p
                for j, a in enumerate(v):
                    if j:
                        for shift, coefficient in pp:
                            if j - 1 + shift < p:
                                w[j - 1 + shift] += j * a * coefficient
                return [a % p for a in w]
            y = [0,1] + [0] * (p-2)
            dp_y = y
            for _ in range(p):
                dp_y = apply(dp_y)
            scalar = dp_y[0]  # P(0)=1.
            for j in range(p):
                v = [int(i == j) for i in range(p)]
                dv = apply(v)
                dpv = v
                for _ in range(p):
                    dpv = apply(dpv)
                assert dpv == [(scalar*a) % p for a in dv], (p,power_exponent,j)
            print('BASE_OPERATOR p='+str(p), 'P=(1+y^2)^'+str(power_exponent), 'c='+str(scalar))


def periods():
    count = 0
    for k in range(1,61):
        onset = 1 + max([p**(2*r-1) for p,r in factors(k)] or [0])
        shift = phi(k)**2
        end = onset + 399 + shift
        seqs = {'A281180': modular_s(k,end)}
        for name,(q,e) in MEMBERS.items():
            seqs[name] = modular_c(k,end,q,e)
        for name,seq in seqs.items():
            for i,want in enumerate(PREFIX[name]):
                n = i+1 if name=='A281180' else i
                assert seq[n] == want % k, (name,k,n,'modular-prefix')
            for n in range(onset,onset+400):
                assert seq[n+shift] == seq[n], (name,k,n,shift,seq[n],seq[n+shift])
                count += 1
        print('PERIOD k='+str(k), 'shift='+str(shift), 'onset='+str(onset),
              '400 consecutive comparisons/member; 5 members OK; generated through n='+str(end))
    assert count == 120000
    print('PERIOD total exact comparisons:',count)


def main():
    fidelity()
    base_operator_check()
    periods()
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
