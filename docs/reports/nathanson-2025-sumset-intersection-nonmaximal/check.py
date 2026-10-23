#!/usr/bin/env python3
"""Checks for M. B. Nathanson, "Intersections of sumsets in additive number theory", arXiv:2512.23574v3, Problem 2:
for A a nonbasis of order h but not a maximal nonbasis, is there an asymptotically strictly decreasing (A_q) with
A = cap A_q and hA = cap hA_q?  Here hA is the set of sums of exactly h elements of A (repetition allowed).

Window checks on Z (integers in [-W, W], sums restricted to a smaller window where they are complete).
Fidelity: the introduction's example (A#_q = {|r| >= q} has hA#_q = Z for h >= 2), and the statement that h*Z is a maximal
nonbasis of order h exactly when h is prime (residue computation, h <= 40).
Claim: A_F = 2Z minus a finite nonempty set F of even integers is a nonbasis of order h (hA_F = 2Z) and not maximal; for
sample admissible sequences A_q = A_F with extra odd elements (and part of F), hA_q contains every odd integer of the window,
so cap hA_q contains 1, which is not in hA_F. Exits 1 on failure.
"""
import itertools, random, sys
ok = True
W = 400; INNER = 60
def hsum(A, h, lo=-INNER, hi=INNER):
    A = sorted(A); S = set(A)
    cur = set(A)
    for _ in range(h - 1):
        cur = {x + a for x in cur for a in A if -3 * W <= x + a <= 3 * W}
    return {x for x in cur if lo <= x <= hi}
win = set(range(-INNER, INNER + 1))
# fidelity 1: A#_q
fid1 = all(hsum({r for r in range(-W, W + 1) if abs(r) >= q}, h) == win for q in (1, 3, 7) for h in (2, 3, 4))
print(f"fidelity: A#_q = {{|r| >= q}} has hA#_q = Z on the window for q in (1, 3, 7), h in (2, 3, 4): {fid1}")
# fidelity 2: h*Z maximal nonbasis iff h prime (residues: h-fold sums of hZ u {x} hit residues k*x mod h, k = 0..h)
def is_prime(n): return n > 1 and all(n % p for p in range(2, int(n ** 0.5) + 1))
def maximal_hZ(h):
    return all(len({(k * x) % h for k in range(h + 1)}) == h for x in range(1, h))
fid2 = all(maximal_hZ(h) == is_prime(h) for h in range(2, 41))
print(f"fidelity: h*Z is a maximal nonbasis of order h iff h is prime (h = 2..40, residue form): {fid2}")
ok &= fid1 and fid2
evens = set(range(-INNER, INNER + 1, 2)); odds = win - evens
rng = random.Random(20261024)
claim_ok = True
for F in ({2}, {0}, {2, 4, 10}, {-6, 8}):
    AF = {r for r in range(-W, W + 1) if r % 2 == 0 and r not in F}
    for h in (2, 3, 4):
        nb = hsum(AF, h) == evens
        notmax = hsum(set(range(-W, W + 1, 2)), h) == evens  # the proper superset 2Z is still a nonbasis
        claim_ok &= nb and notmax
        for trial in range(3):
            extra_odd = {r for r in range(-W, W + 1) if r % 2 and rng.random() < 0.05 and abs(r) > 2 * INNER} | {W - 1}
            Fpart = set(list(F)[: rng.randint(0, len(F))])
            Aq = AF | extra_odd | Fpart
            hs = hsum(Aq, h)
            claim_ok &= odds <= hs and 1 in hs and 1 not in hsum(AF, h)
print(f"claim: for F in ({{2}}, {{0}}, {{2,4,10}}, {{-6,8}}) and h in (2,3,4), hA_F = 2Z, A_F is not maximal, and every sample admissible stage A_q has every odd integer of the window in hA_q (so 1 in cap hA_q, 1 not in hA_F): {claim_ok}")
ok &= claim_ok
print("ALL_CHECKS_OK" if ok else "CHECK_FAIL")
sys.exit(0 if ok else 1)
