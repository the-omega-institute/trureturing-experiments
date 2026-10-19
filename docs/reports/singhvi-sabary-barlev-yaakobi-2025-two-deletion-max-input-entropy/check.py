#!/usr/bin/env python3
"""Singhvi-Sabary-Bar-Lev-Yaakobi (arXiv:2407.10026v2, IEEE TIT 71 (2025)), Section VII-A: maximal input entropy
of the binary 2-deletion and 2-insertion channels under uniform input.

Phi_del(y) = sum_x f(N(x, y)) over inputs x of length |y|+2, Phi_ins(y) = sum_x f(N(y, x)) over inputs x of length
|y|-2, f(t) = t log2 t, N = embedding count; maximising the input entropy is minimising Phi.

(1) Exhaustive maximisers: 2-deletion m = 1..12 (only the two alternating words); 2-insertion m = 2..16 (all words tie
    at m = 2, the nonconstant words at m = 3, exactly 0110 and 1001 at m = 4, only the alternating words for m >= 5).
(2) Closed forms Phi_del(A_m) = 4m^2 + 3m log2 3 + (m+1) log2(m+1) and Phi_ins(A_m) = (m-1) log2(m-1).
(3) Run formulas: the multiset of 2-deletion embedding counts (adjacent-run components) and the run expression for
    Phi_del, checked against brute force for every word of length <= 12.
(4) The 2-deletion minimum is attained only by the constant words (the source's minimisation theorem), m = 1..12.
Exit 0 and RESULT PASS when every check holds.
"""
import itertools, math, sys
from collections import Counter

def f(t): return t * math.log2(t) if t > 0 else 0.0
def runs(y): return [len(list(g)) for _, g in itertools.groupby(y)]
def hist_del(y):
    n = len(y) + 2; c = Counter()
    for i in range(n):
        for j in range(i + 1, n):
            for a in '01':
                for b in '01':
                    c[y[:i] + a + y[i:j - 1] + b + y[j - 1:]] += 1
    return c
def hist_ins(y):
    m = len(y); c = Counter()
    for i in range(m):
        for j in range(i + 1, m): c[y[:i] + y[i + 1:j] + y[j + 1:]] += 1
    return c
def phi(c): return sum(f(v) for v in c.values())
def words(m): return [format(t, '0%db' % m) for t in range(2 ** m)] if m else ['']
def alt(m): return ('01' * m)[:m]
def comp(y): return ''.join('1' if ch == '0' else '0' for ch in y)
def argmin(vals): mn = min(vals.values()); return sorted(y for y, v in vals.items() if abs(v - mn) < 1e-9)
def ins_multiset(y):
    r = runs(y); R = len(r); out = [math.comb(x, 2) for x in r if x >= 2]
    out += [r[i] * r[j] for i in range(R) for j in range(i + 2, R)]
    comps = []; cur = None
    for e in range(R - 1):
        w = r[e] * r[e + 1]
        if cur is None: cur = w
        elif r[e] == 1: cur += w
        else: comps.append(cur); cur = w
    if cur is not None: comps.append(cur)
    return sorted(out + comps)
def phi_del_runs(y):
    m = len(y); r = runs(y)
    s = sum(f(math.comb(x + 2, 2)) + 2 * (m + 1 - x) * f(x + 1) for x in r)
    longs = [k for k, x in enumerate(r) if x >= 2]
    if not longs: return s + f(m + 1)
    u = [r[k] for k in longs]
    gaps = [longs[t + 1] - longs[t] - 1 for t in range(len(u) - 1)]
    C = sum((x - 2) * f(x + 1) for x in u) + f(u[0] + longs[0] + 1) + f(u[-1] + len(r) - longs[-1])
    C += sum(f(u[t] + u[t + 1] + gaps[t] + 1) for t in range(len(u) - 1))
    return s + C

ok = True
def check(cond, msg):
    global ok
    print(('PASS ' if cond else 'FAIL ') + msg)
    ok = ok and cond
for m in range(1, 13):
    vals = {y: phi(hist_del(y)) for y in words(m)}
    check(argmin(vals) == sorted([alt(m), comp(alt(m))]), f'deletion m={m}: maximisers are the alternating words')
    check(abs(vals[alt(m)] - (4 * m * m + 3 * m * math.log2(3) + f(m + 1))) < 1e-6, f'deletion m={m}: closed form')
    mx = max(vals.values())
    check(sorted(y for y, v in vals.items() if abs(v - mx) < 1e-9) == sorted(['0' * m, '1' * m]), f'deletion m={m}: minimisers are constant')
    check(all(abs(vals[y] - phi_del_runs(y)) < 1e-6 for y in vals), f'deletion m={m}: run formula')
    if m >= 2: check(all(sorted(hist_ins(y).values()) == ins_multiset(y) for y in vals), f'insertion m={m}: embedding-count multiset')
for m in range(2, 17):
    vals = {y: phi(hist_ins(y)) for y in words(m)}
    a = argmin(vals)
    if m == 2: exp = sorted(words(2))
    elif m == 3: exp = sorted(y for y in words(3) if len(set(y)) == 2)
    elif m == 4: exp = ['0110', '1001']
    else: exp = sorted([alt(m), comp(alt(m))])
    check(a == exp, f'insertion m={m}: maximisers {a[:6]}')
    if m >= 5: check(abs(vals[alt(m)] - f(m - 1)) < 1e-6, f'insertion m={m}: closed form')
print('RESULT PASS' if ok else 'RESULT FAIL')
sys.exit(0 if ok else 1)
