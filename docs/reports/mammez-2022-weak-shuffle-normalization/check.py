"""Check of the proposed proof of Mammez, Weak stuffle algebras (arXiv:2003.09142v3, IEJA 32 (2022)), Conjecture 3.10.

Weak shuffle product on words: (au) * (bv) = f1(a,b) a(u * bv) + f2(a,b) b(au * v), eps the unit.
Random coefficient systems over alphabets of size 2..4 with values in {0, 1, 2, -3, 1/2} are kept when the
product is commutative and associative on all words up to a length bound. For each kept system with an
exceptional pair (a != b, f1(a,b) not in {0,1}), the normalized system g replaces every such value by 1,
and Phi(w) = q(w) w with q(w) = prod_e f(e)^(-h_e(w)) is checked to be an algebra isomorphism on short words.
"""
import itertools, random
from fractions import Fraction as Fr
from functools import lru_cache
random.seed(310)

def make_product(F1, F2):
    @lru_cache(maxsize=None)
    def prod(u, v):
        if not u: return {v: Fr(1)}
        if not v: return {u: Fr(1)}
        a, b = u[0], v[0]; out = {}
        c1 = F1[(a, b)]
        if c1:
            for w, c in prod(u[1:], v).items():
                out[(a,) + w] = out.get((a,) + w, 0) + c1 * c
        c2 = F2[(a, b)]
        if c2:
            for w, c in prod(u, v[1:]).items():
                out[(b,) + w] = out.get((b,) + w, 0) + c2 * c
        return {w: c for w, c in out.items() if c}
    return prod

def lin(prod, x, y):
    out = {}
    for u, cu in x.items():
        for v, cv in y.items():
            for w, c in prod(u, v).items():
                out[w] = out.get(w, 0) + cu * cv * c
    return {w: c for w, c in out.items() if c}

def words(X, L):
    for l in range(L + 1):
        yield from itertools.product(X, repeat=l)

def is_comm_assoc(prod, X, L):
    W = list(words(X, L))
    for u in W:
        for v in W:
            if prod(u, v) != prod(v, u): return False
    for u in W:
        for v in W:
            for w in W:
                if lin(prod, prod(u, v), {w: Fr(1)}) != lin(prod, {u: Fr(1)}, prod(v, w)): return False
    return True

def h(e, w, F1):
    a, b = e
    P = {x for x in set(w) | {a, b} if F1.get((x, b), 0) != 0 and x != b}
    # P is defined over the whole alphabet: letters x with f(x, b) != 0
    i = 0
    while i < len(w) and w[i] in PP[e]: i += 1
    seg = w[:i]; nxt = w[i] if i < len(w) else None
    return seg.count(a) if nxt == b else 0

VALS = [Fr(0), Fr(1), Fr(2), Fr(-3), Fr(1, 2)]
kept = exc = 0; bad = 0; byn = {}
for trial in range(30000):
    n = random.choice([2, 3, 3, 4])
    X = tuple(range(n))
    F1 = {}; F2 = {}
    for a in X:
        d = random.choice([Fr(0), Fr(1)])
        F1[(a, a)] = d; F2[(a, a)] = d
    for a in X:
        for b in X:
            if a != b:
                F1[(a, b)] = random.choices([Fr(0), Fr(1), random.choice([Fr(2), Fr(-3), Fr(1, 2)])], weights=[4, 4, 3])[0]
    for a in X:
        for b in X:
            if a != b: F2[(a, b)] = F1[(b, a)]
    prod = make_product(F1, F2)
    if not is_comm_assoc(prod, X, 2): continue
    kept += 1
    E = [(a, b) for a in X for b in X if a != b and F1[(a, b)] not in (0, 1)]
    if not E: continue
    exc += 1; byn[n] = byn.get(n, 0) + 1
    G1 = {k: (Fr(1) if (k[0] != k[1] and v not in (0, 1)) else v) for k, v in F1.items()}
    G2 = {(a, b): (G1[(b, a)] if a != b else F2[(a, b)]) for a in X for b in X}
    gprod = make_product(G1, G2)
    PP = {e: {x for x in X if F1[(x, e[1])] != 0 and x != e[1]} for e in E}
    def q(w):
        r = Fr(1)
        for e in E:
            r *= F1[e] ** (-h(e, w, F1))
        return r
    W = list(words(X, 3))
    for u in W:
        for v in W:
            lhs = {w: q(w) * c for w, c in prod(u, v).items()}
            rhs = {w: q(u) * q(v) * c for w, c in gprod(u, v).items()}
            if lhs != rhs:
                bad += 1; print("FAIL", n, E, u, v); break
        else: continue
        break
    if not is_comm_assoc(gprod, X, 2): bad += 1; print("NORMALIZED NOT ASSOC", E)
print("systems kept (comm+assoc on words <= 2):", kept, "with exceptional pairs:", exc, "by alphabet size:", dict(sorted(byn.items())))
if exc < 20 or byn.get(4, 0) < 20: bad += 1; print("FAIL too few exceptional systems")
print("bad=", bad)
raise SystemExit(1 if bad else 0)
