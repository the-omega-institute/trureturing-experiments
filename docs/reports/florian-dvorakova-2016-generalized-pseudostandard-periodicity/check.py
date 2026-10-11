#!/usr/bin/env python3
"""Florian-Dvorakova (EJC 23(1) P1.2, 2016) Conjecture 29: periodicity of d-ary generalized
pseudostandard words. Finite-prefix sanity checks of the proposed proof route (not a proof).

An involutory antimorphism is theta_pi(w) = reverse(map pi w) with pi an involution of the alphabet.
Checks:
  A. Condition 2 of the conjecture is equivalent to pi_{i+2} = pi_{i+1} pi_i pi_{i+1}
     (exhaustive over alphabets of size 2..4, words up to length 5).
  B. Normalization: between consecutive pseudopalindromic prefixes U[:m] < U[:m'] of a generated word,
     some involution pi gives closure_pi(U[:m] + U[m]) == U[:m'] (the insertion lemma), on random
     directive bi-sequences with periodic tails.
  C. For d-ary words whose normalized directive is observed over enough steps: purely periodic with a
     small period <=> both conditions hold on a tail; and, when they hold, consecutive gaps are
     nonincreasing and the least period equals g * ord(pi_{i+1} pi_i).
Exit 0 and final line `bad= 0` when every check agrees.
"""
import itertools, random, sys

def involutions(d):
    out = []
    for p in itertools.permutations(range(d)):
        if all(p[p[a]] == a for a in range(d)):
            out.append(p)
    return out

def theta(pi, w):
    return tuple(pi[a] for a in reversed(w))

def closure(pi, w):
    # w = p s with s the longest theta-palindromic suffix; closure = p s theta(p)
    for k in range(len(w) + 1):
        s = w[k:]
        if theta(pi, s) == s:
            return w + theta(pi, w[:k])
    raise AssertionError

def closure_brute(pi, w):
    # shortest theta-palindrome with prefix w, by search over extensions of length <= |w|
    d = len(pi)
    for ext in range(len(w) + 1):
        for tail in itertools.product(range(d), repeat=ext):
            z = w + tail
            if theta(pi, z) == z:
                return z
    raise AssertionError

def compose(p, q):  # (p o q)(a) = p(q(a))
    return tuple(p[q[a]] for a in range(len(p)))

def order(p):
    k, cur, ident = 1, p, tuple(range(len(p)))
    while cur != ident:
        cur, k = compose(p, cur), k + 1
    return k

bad = 0
rng = random.Random(20261012)

# A. condition 2 <=> permutation identity
countA = 0
for d in range(2, 5):
    inv = involutions(d)
    pals = {}
    for pi in inv:
        pals[pi] = [w for L in range(0, 6) for w in itertools.product(range(d), repeat=L) if theta(pi, w) == w]
    for p0, p1, p2 in itertools.product(inv, repeat=3):
        literal = all(theta(p2, theta(p1, w)) == theta(p1, w) for w in pals[p0])
        identity = p2 == compose(compose(p1, p0), p1)
        countA += 1
        if literal != identity:
            bad += 1; print('A mismatch', d, p0, p1, p2)
print('A triples', countA)

# closure formula agrees with brute-force shortest closure
for _ in range(400):
    d = rng.randint(2, 4); pi = rng.choice(involutions(d))
    w = tuple(rng.randrange(d) for _ in range(rng.randint(0, 6)))
    if closure(pi, w) != closure_brute(pi, w):
        bad += 1; print('closure mismatch', pi, w)

def generate(delta, Theta, maxlen):
    W, prefs = (), [()]
    n = 0
    while len(W) < maxlen:
        W = closure(Theta(n), W + (delta(n),))
        prefs.append(W); n += 1
        if n > 100000: break
    return W, prefs

def is_pseudopal(w, inv):
    return [pi for pi in inv if theta(pi, w) == w]

def least_period(U, cap):
    for p in range(1, cap + 1):
        if all(U[i + p] == U[i] for i in range(len(U) - p)):
            return p
    return None

statsB = statsC = periodic_cases = aperiodic_cases = skipped = long_period_cases = 0
for trial in range(1500):
    d = rng.randint(2, 5); inv = involutions(d)
    r, q = rng.randint(0, 3), rng.randint(1, 4)
    pre = [(rng.randrange(d), rng.choice(inv)) for _ in range(r)]
    per = [(rng.randrange(d), rng.choice(inv)) for _ in range(q)]
    seq = lambda n: pre[n] if n < r else per[(n - r) % q]
    L = 1500
    U, _ = generate(lambda n: seq(n)[0], lambda n: seq(n)[1], L)
    U = U[:L]
    if len(set(U)) < d:
        skipped += 1; continue
    lens = [m for m in range(0, L + 1) if is_pseudopal(U[:m], inv)]
    # B: insertion lemma between consecutive pseudopalindromic prefixes
    pis = []
    for a, b in zip(lens, lens[1:]):
        cands = [pi for pi in inv if theta(pi, U[:b]) == U[:b] and closure(pi, U[:a] + (U[a],)) == U[:b]]
        statsB += 1
        if not cands:
            bad += 1; print('B no closure step', d, a, b); break
        pis.append(cands)
    else:
        # C: periodicity vs tail conditions on the normalized directive
        steps = len(pis)
        if steps < 12:
            continue
        p = least_period(U, L // 8)
        full = [i for i in range(steps) if len(set(U[:lens[i + 1]])) == d]
        if len(full) < 10:
            continue
        uniq = [pis[i][0] if len(pis[i]) == 1 else None for i in range(steps)]
        tail = full[len(full) // 2:]
        delta_n = [U[lens[i]] for i in range(steps)]  # letter appended at step i (to reach lens[i+1])
        def conds(ix):
            # i in ix: pi_i = uniq[i] fixes U[:lens[i+1]]; next letter delta_{i+1} = U[lens[i+1]]
            ok1 = len({uniq[i][delta_n[i + 1]] for i in ix if i + 1 < steps}) == 1
            ok2 = all(uniq[i + 2] == compose(compose(uniq[i + 1], uniq[i]), uniq[i + 1]) for i in ix if i + 2 < steps)
            return ok1 and ok2
        if any(uniq[i] is None for i in tail):
            bad += 1; print('C nonunique antimorphism on all-letter prefix'); continue
        c = conds(tail)
        statsC += 1
        if p is not None:
            periodic_cases += 1
            if not c:
                bad += 1; print('C periodic but conditions fail', d, p)
            else:
                gaps = [lens[i + 2] - lens[i + 1] for i in tail if i + 2 < len(lens)]
                if any(x < y for x, y in zip(gaps, gaps[1:])):
                    bad += 1; print('C gaps increase', gaps[:10])
                g = gaps[-1]
                rho = compose(uniq[tail[-2]], uniq[tail[-3]])
                if g * order(rho) != p:
                    bad += 1; print('C period formula', g, order(rho), p)
        else:
            if c:
                # conditions hold: the predicted period g * ord(rho) may exceed the small-period cap
                gaps = [lens[i + 2] - lens[i + 1] for i in tail if i + 2 < len(lens)]
                g = gaps[-1]
                rho = compose(uniq[tail[-2]], uniq[tail[-3]])
                P = g * order(rho)
                if P < len(U) and all(U[i + P] == U[i] for i in range(len(U) - P)) and least_period(U, P) == P:
                    long_period_cases += 1
                else:
                    bad += 1; print('C conditions hold but no period g*ord(rho)', d, g, order(rho), len(tail))
            else:
                aperiodic_cases += 1

print('B steps', statsB, 'C cases', statsC, 'periodic', periodic_cases, 'long-period periodic', long_period_cases, 'aperiodic', aperiodic_cases, 'not d-ary', skipped)
print('bad=', bad)
sys.exit(0 if bad == 0 else 1)
