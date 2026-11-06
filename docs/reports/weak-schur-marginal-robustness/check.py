#!/usr/bin/env python3
"""Exact checks for fixed-k weak-Schur marginal robustness.

All evidential arithmetic uses integers or Fraction. RSK shape probabilities on
classical iid words equal diagonal-input WSS probabilities. The independent
n=2,3 character-projector calculation checks this identification in small cases.
No finite experiment asserts the n-to-infinity conclusion.
"""
from collections import Counter, defaultdict
from fractions import Fraction as Q
from itertools import combinations, product
from math import comb, factorial, prod


def shape(word):
    rows = []
    for letter in word:
        for row in rows:
            j = next((j for j, x in enumerate(row) if x > letter), len(row))
            if j == len(row):
                row.append(letter)
                break
            row[j], letter = letter, row[j]
        else:
            rows.append([letter])
    return tuple(map(len, rows))


def law(weights, n):
    out = defaultdict(int)
    for word in product(range(len(weights)), repeat=n):
        out[shape(word)] += prod(weights[i] for i in word)
    den = sum(weights) ** n
    return {s: Q(w, den) for s, w in out.items()}


def tv(p, q):
    return sum((abs(a-b) for a, b in zip(p, q)), Q(0)) / 2


def loss(s, n, target):
    padded = s + (0,) * (len(target)-len(s))
    return tv([Q(a, n) for a in padded], target)


def integer_partitions(n, lo=1):
    if n == 0:
        yield ()
    for first in range(lo, n+1):
        for rest in integer_partitions(n-first, first):
            yield (first,) + rest


def word_probability(counts, d, r):
    # A fixed word has m distinct symbols; B records which are in S.
    m = len(counts)
    numerator = 0
    for mask in range(1 << m):
        j = mask.bit_count()
        if 0 <= r-j <= d-m:
            exponent = sum(c for i, c in enumerate(counts) if mask >> i & 1)
            numerator += comb(d-m, r-j) * 3**exponent
    return Q(numerator, comb(d, r) * (2*d)**sum(counts))


def multiplicity(counts, d):
    m = len(counts)
    falling = factorial(d) // factorial(d-m)
    return (falling * factorial(sum(counts)) //
            (prod(factorial(c) for c in counts) *
             prod(factorial(v) for v in Counter(counts).values())))


def marginal_tv(d, k):
    r = d//2
    total, mass, words = Q(0), Q(0), 0
    for counts in integer_partitions(k):
        if len(counts) > d:
            continue
        mult = multiplicity(counts, d)
        p = word_probability(counts, d, r)
        mass += mult*p
        words += mult
        total += mult*abs(p-Q(1, d**k)) / 2
    assert mass == 1 and words == d**k
    return total


def run():
    # Fidelity: source Proposition “Marginal constraints admit exact failure
    # at epsilon=0”, parts (i),(ii): tau=sum_i p_i |i><i|^tensor n.
    for p in [(Q(3, 4), Q(1, 4)), (Q(1, 2), Q(1, 3), Q(1, 6)),
              (Q(1, 4),)*4]:
        for n in range(1, 9):
            L = defaultdict(Q)
            marginal = [Q(0)]*len(p)
            for i, mass in enumerate(p):
                s = shape((i,)*n)
                L[s] += mass
                marginal[i] += mass
            assert dict(L) == {(n,): Q(1)}
            assert marginal == list(p)
            assert sum(a*loss(s, n, p) for s, a in L.items()) == 1-p[0]
    print('FIDELITY: source marginal-failure proposition (i),(ii), 24 exact cases OK')

    # Independent central-character formula: P_lambda = f_lambda/n! sum_pi
    # chi_lambda(pi) U_pi, with traces of cycles = power sums of the spectrum.
    for weights in [(3, 3, 1, 1), (1, 1, 1, 1), (3, 1), (4, 2, 1)]:
        p = [Q(w, sum(weights)) for w in weights]
        t2, t3 = sum(x*x for x in p), sum(x**3 for x in p)
        l2 = law(weights, 2)
        assert l2.get((2,), 0) == (1+t2)/2
        assert l2.get((1, 1), 0) == (1-t2)/2
        l3 = law(weights, 3)
        assert l3.get((3,), 0) == (1+3*t2+2*t3)/6
        assert l3.get((2, 1), 0) == 2*(1-t3)/3
        assert l3.get((1, 1, 1), 0) == (1-3*t2+2*t3)/6
    print('WSS CROSS-CHECK: n=2,3 character projectors agree with exact RSK laws')

    cases = 0
    for r in range(2, 17):
        d = 2*r
        for k in range(2, min(r, 8)+1):
            defect = marginal_tv(d, k)
            assert defect <= Q(k*(k-1), 4*d), (d, k, defect)
            cases += 1
    print(f'MARGINAL BOUND: k(k-1)/(4d), {cases} exact cases, 4<=d<=32, 2<=k<=min(d/2,8) OK')
    for r in range(2, 5):
        d = 2*r
        for k in range(1, min(r, 4)+1):
            brute = sum((abs(word_probability(tuple(Counter(w).values()), d, r)-Q(1,d**k))
                         for w in product(range(d), repeat=k)), Q(0))/2
            assert brute == marginal_tv(d, k)
    print('MARGINAL CROSS-CHECK: direct word enumeration agrees with occupancy multiplicities')

    d, r = 4, 2
    sets = list(combinations(range(d), r))
    assert len(sets) == 6
    pair = {}
    for x, y in product(range(d), repeat=2):
        pair[x, y] = sum((Q(3 if x in S else 1, 8)*Q(3 if y in S else 1, 8)
                         for S in sets), Q(0))/6
        assert pair[x, y] == (Q(5,64) if x == y else Q(11,192))
    assert sum(pair.values()) == 1
    assert sum((abs(x-Q(1,16)) for x in pair.values()),Q(0))/2 == Q(1,16)
    print('EXPLICIT d=4: 6 components; pair diagonal=5/64, off-diagonal=11/192; TD=1/16')

    for r in range(2, 17):
        d = 2*r
        spectrum = [Q(3,2*d)]*r + [Q(1,2*d)]*r
        assert min(spectrum) > 0 and sum(spectrum) == 1
        assert tv(spectrum, [Q(1,d)]*d) == Q(1,4)
        for n in range(1, 9):
            assert min(spectrum)**n == Q(1,(2*d)**n) > 0
    print('SPECTRUM/RANK: 3/(2d),1/(2d) each repeated d/2; TV=1/4; positive tensor lower bound OK')
    for d in (4, 6):
        r = d//2
        for n in range(1, 6):
            expected = law((3,)*r+(1,)*r, n)
            mixture = defaultdict(Q)
            for S in combinations(range(d), r):
                component = law(tuple(3 if i in S else 1 for i in range(d)), n)
                assert component == expected
                for s, mass in component.items():
                    mixture[s] += mass/Q(comb(d,r))
            assert dict(mixture) == expected and sum(mixture.values()) == 1
    print('WSS MIXTURE: all 6 (d=4) and 20 (d=6) components and their mean have identical laws, 1<=n<=5')
    for n in range(1, 8):
        L = law((3,3,1,1), n)
        err = sum((a*loss(s,n,[Q(1,4)]*4) for s,a in L.items()),Q(0))
        own = sum((a*loss(s,n,[Q(3,8),Q(3,8),Q(1,8),Q(1,8)]) for s,a in L.items()),Q(0))
        assert abs(err-Q(1,4)) <= own
        print(f'd=4 n={n}: target_error={err}; own_spectrum_error={own}; reverse_triangle_OK')
    print('LIMIT SCOPE: finite checks only; fixed-d noiseless WSS consistency is a required proof ingredient')
    print('ALL_CHECKS_OK')

if __name__ == '__main__':
    run()
