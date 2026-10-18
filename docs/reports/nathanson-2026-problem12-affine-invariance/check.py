"""Nathanson, arXiv:2604.04781, Problem 12: if A'_q = sA_q + t (pointwise, s not in {0, 1, -1}), is H(A_q) = H(A'_q)?
Here H(A_q) = {h >= 1 : hA = intersection of hA_q}, A = intersection of A_q, and hB is the h-fold sumset.

The proposed answer is yes, for every nonzero s. The map x -> s x + h t is injective, it commutes with intersections,
and h(sB + t) = s(hB) + ht. This program checks the identity on random finite families of finite integer sets (not necessarily nested,
since the identity holds for every family), for s in {-3, -2, 2, 3, 5}, several t and 1 <= h <= 4. It also reports, as a control, how often H misses some h,
so that the check is not vacuous. Standard library only.
"""
import random, itertools
random.seed(20261019)

def sumset(B, h):
    S = {0}
    for _ in range(h):
        S = {x + b for x in S for b in B}
    return S

def H(fam, hmax):
    A = set.intersection(*fam)
    return {h for h in range(1, hmax + 1) if sumset(A, h) == set.intersection(*(sumset(B, h) for B in fam))}

def random_family():
    # Arbitrary (not nested) finite families sharing a common core: the identity is checked for every family,
    # and non-nested families make H(A_q) nontrivial.
    core = set(random.sample(range(-6, 7), random.randint(1, 3)))
    return [core | set(random.sample(range(-12, 13), random.randint(1, 6))) for _ in range(random.randint(2, 5))]

cases = mismatches = nontrivial = 0
for _ in range(400):
    fam = random_family()
    if not set.intersection(*fam):
        continue
    h0 = H(fam, 4)
    nontrivial += h0 != set(range(1, 5))
    for s in (-3, -2, 2, 3, 5):
        for t in (-4, 0, 7):
            fam2 = [{s * a + t for a in B} for B in fam]
            cases += 1
            mismatches += H(fam2, 4) != h0
print("cases:", cases, "mismatches:", mismatches, "families where H misses some h <= 4:", nontrivial)
ok = cases > 0 and mismatches == 0 and nontrivial > 0
print("RESULT", "PASS" if ok else "FAIL")
raise SystemExit(0 if ok else 1)
