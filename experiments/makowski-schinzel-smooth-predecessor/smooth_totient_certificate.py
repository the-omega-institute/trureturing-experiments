#!/usr/bin/env python3
"""Exact certificate: sigma(phi(n)) >= n/2 for every n all of whose prime factors q have B-smooth q-1.

Reduction: it suffices that every odd squarefree m with this property satisfies
    prod_{q | m, q not Fermat} q/(q-1) <= prod_{odd p <= B} s_p(e_p),   s_p(e) = sigma(p^e)/p^e,
since the Fermat part s_2(A) prod_{Fermat q | m}(1 - 1/q) is at least 1.

Certificate (Cohen-type grouping with a complete table and an Euler-product tail):
  * candidates: every q = m + 1 with m <= X even, B-smooth, not a power of two, and not proved
    composite (Miller-Rabin in gmpy2 never declares a prime composite; survivors need no proof);
  * each candidate gets a label g(q), an odd prime dividing q - 1; candidates above X are labelled
    by the largest odd prime factor of q - 1;
  * for each label p, with members c_1 < c_2 < ... and L_j = prod_{i<=j} (1 - 1/c_i), check
        s_p(j) L_j >= 1  (j >= 1)   and   s_p(j+1) L_j (1 - T_p) >= 1  (j >= 0),
    where T_p = (1/(2p)) (prod_{r<=p} r/(r-1) - sum_{m' <= X/(2p), m' p-smooth} 1/m') bounds the
    reciprocal mass of labelled primes above X.
All comparisons use integers: L_j is a floor-rounded fixed-point lower bound with scale 10^60, and the
finite part of T_p is a floor-rounded lower bound, so T_p is over-estimated.

Usage: python3 smooth_totient_certificate.py B X [--table FILE]
"""
import argparse, json, math, sys
from fractions import Fraction
import gmpy2
from sympy import primerange

SCALE = 10 ** 60


def smooth(primes, limit):
    out = [1]
    for p in primes:
        new = []
        for m in out:
            x = m
            while x <= limit:
                new.append(x)
                x *= p
        out = new
    return sorted(out)


def log_s(p, k):
    return math.log((1 - p ** -(k + 1)) / (1 - 1 / p))


def margin(p, members, log_one_minus_tail):
    acc, best = 0.0, math.inf
    for j in range(len(members) + 1):
        if j > 0:
            best = min(best, log_s(p, j) + acc)
        best = min(best, log_s(p, j + 1) + acc + log_one_minus_tail)
        if j < len(members):
            acc += math.log1p(-1 / members[j])
    return best


def tails(primes, odd, limit):
    out = {}
    for p in odd:
        euler = Fraction(1)
        for r in primes:
            if r <= p:
                euler *= Fraction(r, r - 1)
        finite = sum(SCALE // (2 * p * m) for m in smooth([r for r in primes if r <= p], limit // (2 * p)))
        out[p] = euler / (2 * p) - Fraction(finite, SCALE)
    return out


def assign(cands, odd, tail, iters=20000):
    lt = {p: math.log1p(-float(tail[p])) for p in odd}
    cap = {p: math.log(p / (p - 1)) + lt[p] for p in odd}
    load = {p: 0.0 for p in odd}
    groups = {p: [] for p in odd}
    allowed = {q: [p for p in odd if (q - 1) % p == 0] for q in cands}
    order = sorted(q for q in cands if len(allowed[q]) == 1) + sorted(q for q in cands if len(allowed[q]) > 1)
    for q in order:
        p = max(allowed[q], key=lambda r: cap[r] - load[r])
        load[p] += -math.log1p(-1 / q)
        groups[p].append(q)
    for p in odd:
        groups[p].sort()
    mg = {p: margin(p, groups[p], lt[p]) for p in odd}
    for _ in range(iters):
        g = min(odd, key=lambda p: mg[p])
        if mg[g] > 2e-9:
            break
        best = None
        for q in groups[g]:
            w = -math.log1p(-1 / q)
            for h in allowed[q]:
                if h != g:
                    est = min(mg[g] + w, mg[h] - w)
                    if best is None or est > best[0]:
                        best = (est, q, h)
        if best is None or best[0] <= mg[g]:
            break
        _, q, h = best
        groups[g].remove(q)
        groups[h].append(q)
        groups[h].sort()
        mg[g] = margin(g, groups[g], lt[g])
        mg[h] = margin(h, groups[h], lt[h])
    return groups


def check(groups, odd, tail):
    worst = None
    for p in odd:
        fixed = SCALE
        one_minus_tail = 1 - tail[p]
        members = groups[p]
        for j in range(len(members) + 1):
            for k, extra in ((j, Fraction(1)), (j + 1, one_minus_tail)):
                if k == 0:
                    continue
                ratio = Fraction(p ** (k + 1) - 1, (p - 1) * p ** k) * fixed * extra / SCALE
                if worst is None or ratio < worst[0]:
                    worst = (ratio, p, j, k)
            if j < len(members):
                q = members[j]
                fixed = fixed * (q - 1) // q
    return worst


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("B", type=int)
    ap.add_argument("X", type=float)
    ap.add_argument("--table")
    a = ap.parse_args()
    B, X = a.B, int(a.X)
    primes = list(primerange(2, B + 1))
    odd = primes[1:]
    cands = [m + 1 for m in smooth(primes, X) if m % 2 == 0 and m & (m - 1) and gmpy2.is_prime(m + 1)]
    tail = tails(primes, odd, X)
    assert all(0 <= tail[p] < 1 for p in odd)
    groups = assign(cands, odd, tail)
    worst = check(groups, odd, tail)
    ok = worst[0] >= 1
    print(json.dumps({"B": B, "X": X, "candidates": len(cands), "result": "PASS" if ok else "FAIL",
                      "min_ratio": float(worst[0]), "at": {"p": worst[1], "j": worst[2], "k": worst[3]}}))
    if a.table and ok:
        label = {q: p for p, qs in groups.items() for q in qs}
        with open(a.table, "w") as f:
            json.dump({"B": B, "X": X, "table": [[q, label[q]] for q in sorted(cands)]}, f)
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
