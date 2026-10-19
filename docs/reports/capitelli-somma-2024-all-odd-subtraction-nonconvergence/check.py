#!/usr/bin/env python3
"""Capitelli-Somma (arXiv:2405.19593v1), Section 3: an all-odd finite subtraction set gives a nonconvergent winning
probability sequence a_n^S (a_n = 0 without a legal move, else 1 - average of a_{n-s} over legal moves s <= n).

Exact rational checks over every all-odd S contained in {1, 3, ..., 15} with |S| <= 4 (162 sets), n < max S + 400:
(1) d_n = b_n + b_{n+1} >= 0 and d_{min S - 1} = 2, where b_n = (-1)^n (1 - 2 a_n);
(2) I_n = sum_{s in S} sum_{r < s} d_{n-r} is constant for n >= max S - 1 and at least 2;
(3) the tail oscillation |a_n - a_{n+1}| stays bounded away from 0 (nonconvergence);
(4) alpha_1 = -I/(4 sum S) predicts the even/odd limits 1/2 +- alpha_1 (checked for five gcd-one sets);
(5) the source's closed form for alpha_1 (Corollary, Eq. alpha_1Final, |S| = 2 odd coprime) is compared with the exact value:
    it agrees for k = 1 and disagrees for every listed pair with k >= 3 (reported, not a failure of this check).
Exit 0 and RESULT PASS when (1)-(4) hold.
"""
from fractions import Fraction as F
from math import gcd
import itertools, sys
def seq(S, N):
    S = sorted(S); a = []
    for n in range(N):
        A = [s for s in S if s <= n]
        a.append(F(0) if not A else 1 - sum(a[n - s] for s in A) / len(A))
    return a
def invariant(S, a):
    b = [(-1) ** n * (1 - 2 * x) for n, x in enumerate(a)]; d = [b[n] + b[n + 1] for n in range(len(b) - 1)]
    return b, d
ok = True; nsets = 0
for t in range(1, 5):
    for S in itertools.combinations(range(1, 16, 2), t):
        nsets += 1; M = max(S); m = min(S); N = M + 400
        a = seq(S, N); b, d = invariant(S, a)
        if any(x < 0 for x in d) or d[m - 1] != 2: ok = False; print('(1) fails', S)
        I = {sum(d[n - r] for s in S for r in range(s)) for n in range(M - 1, N - 1)}
        if len(I) != 1 or min(I) < 2: ok = False; print('(2) fails', S)
        if max(abs(a[n] - a[n + 1]) for n in range(N - 50, N - 1)) < F(1, 1000): ok = False; print('(3) fails', S)
print('(1)-(3) checked on', nsets, 'all-odd sets')
for S in [(1, 3), (3, 5), (1, 3, 5), (3, 5, 7), (5, 7, 9, 11)]:
    a = seq(S, 3000); b, d = invariant(S, a); M = max(S)
    al = -sum(d[M - 1 - r] for s in S for r in range(s)) / (4 * sum(S))
    good = abs(float(a[2998] - (F(1, 2) + al))) < 1e-9 and abs(float(a[2999] - (F(1, 2) - al))) < 1e-9
    print('(4)', S, 'alpha_1 =', al, 'even/odd limits predicted:', good); ok = ok and good
def paper(k, l):
    q = l // k
    if q % 2 == 0: return F(1, l + k) * (-F(1, 2) * q - k * (q + 1) + l - 1)
    return F(1, l + k) * (-F(1, 2) * (q - 1) - l + k * q - 1)
agree = disagree = 0; agree_k = set(); disagree_k = set()
for k in range(1, 16, 2):
    for l in range(k + 2, 30, 2):
        if gcd(k, l) != 1: continue
        a = seq((k, l), l + 5); b, d = invariant((k, l), a)
        true = -sum(d[l - 1 - r] for s in (k, l) for r in range(s)) / (4 * (k + l))
        if true == paper(k, l): agree += 1; agree_k.add(k)
        else:
            disagree += 1; disagree_k.add(k)
            if disagree <= 3: print('(5) source Corollary alpha_1', (k, l), 'source', paper(k, l), 'exact', true)
print('(5) source alpha_1 formula agrees on', agree, 'pairs (k in', sorted(agree_k), '), disagrees on', disagree, 'pairs (k in', sorted(disagree_k), ')')
print('RESULT PASS' if ok else 'RESULT FAIL'); sys.exit(0 if ok else 1)
