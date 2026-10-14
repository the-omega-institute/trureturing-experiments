#!/usr/bin/env python3
"""Checks for N. Bonin, "Conditions of positivity on a shadow Markoff tree", J. Integer Seq. 27 (2024), Article 24.8.3;
arXiv:2403.07668v2, Conjecture 1: the initial shadows (alpha, beta, gamma) of the tree rooted at (1+alpha e, 1+beta e, 1+gamma e)
give a positive integer tree exactly for the rational points of the quadrilateral (0:0:1), (1:0:2), (1:1:1), (0:2:1).
The tree is the source's path(alpha, beta, gamma, C): root (1,1,1), mandatory 'r' then 'l', then any word C (Section 5 code),
with the source's rational mutation formulas.
Claimed reading: for nonzero (alpha, beta, gamma) in N^3, every entry on every path is a nonnegative integer
iff alpha + beta <= 2 gamma and 2 alpha <= beta + gamma.
1. Exact (Fractions), all nonzero triples with alpha, beta, gamma <= 9:
   - every word of length <= 11 gives integral entries;
   - inside the region (both inequalities), every shadow on every word of length <= 11 is nonnegative;
   - outside the region, a negative shadow appears along the all-'l' spine when alpha + beta > 2 gamma, and along the all-'r'
     spine when 2 alpha > beta + gamma (searched to 400 steps; the forced-recurrence mechanism of the proof).
2. The two linear inequalities cut out exactly the quadrilateral (vertex check in the chart gamma = 1).
3. The source's own spot checks: path(1, 0.9, 1, 'rlrlrlr') is negative (Section 5.1).
Prints ALL_OK.
"""
from fractions import Fraction as Fr
import itertools
ok = True
def mut(t, i):
    a, al, b, be, c, ga = t
    if i == 'l':
        nb = (a*a + c*c) / b
        return (a, al, c, ga, nb, (-nb*be + 2*a*al + 2*c*ga) / b)
    na = (b*b + c*c) / a
    return (c, ga, b, be, na, (-na*al + 2*b*be + 2*c*ga) / a)
def tree_ok(al, be, ga, depth):
    """Return (all integral, all shadows nonnegative) over every word of length <= depth."""
    root = (Fr(1), Fr(al), Fr(1), Fr(be), Fr(1), Fr(ga))
    t1 = mut(root, 'r'); t2 = mut(t1, 'l')
    nodes = [root, t1, t2]; frontier = [t2]
    for _ in range(depth):
        nxt = []
        for t in frontier:
            for i in 'lr':
                u = mut(t, i); nodes.append(u); nxt.append(u)
        frontier = nxt
    integral = all(x.denominator == 1 for t in nodes for x in t)
    nonneg = all(t[k] >= 0 for t in nodes for k in (1, 3, 5))
    return integral, nonneg
def spine_negative(al, be, ga, d, steps=400):
    t = (Fr(1), Fr(al), Fr(1), Fr(be), Fr(1), Fr(ga)); t = mut(mut(t, 'r'), 'l')
    for n in range(1, steps + 1):
        t = mut(t, d)
        if min(t[k] for k in (1, 3, 5)) < 0: return n
    return None
checked = inside = outside = 0; latest = 0
for al, be, ga in itertools.product(range(10), repeat=3):
    if (al, be, ga) == (0, 0, 0): continue
    checked += 1
    integral, nonneg = tree_ok(al, be, ga, 11)
    if not integral: print('NONINTEGRAL', al, be, ga); ok = False
    u, v = 2*ga - al - be, be + ga - 2*al
    if u >= 0 and v >= 0:
        inside += 1
        if not nonneg: print('INSIDE_NEGATIVE', al, be, ga); ok = False
    else:
        outside += 1
        n = spine_negative(al, be, ga, 'l') if u < 0 else spine_negative(al, be, ga, 'r')
        if n is None: print('OUTSIDE_NO_NEGATIVE', al, be, ga); ok = False
        else: latest = max(latest, n)
print('checked %d nonzero triples in [0,9]^3: inside %d all nonnegative to depth 11; outside %d all negative on the forced spine (latest at step %d)' % (checked, inside, outside, latest))
# quadrilateral = {x>=0, y>=0, x+y<=2, 2x<=y+1} in the chart gamma=1; vertices
V = [(Fr(0), Fr(0)), (Fr(1, 2), Fr(0)), (Fr(1), Fr(1)), (Fr(0), Fr(2))]
H = [lambda x, y: x >= 0, lambda x, y: y >= 0, lambda x, y: x + y <= 2, lambda x, y: 2*x <= y + 1]
for v in V:
    if not all(h(*v) for h in H): print('VERTEX_FAIL', v); ok = False
tight = [sum(1 for h in [lambda x, y: x == 0, lambda x, y: y == 0, lambda x, y: x + y == 2, lambda x, y: 2*x == y + 1] if h(*v)) for v in V]
if tight != [2, 2, 2, 2]: print('VERTEX_TIGHT_FAIL', tight); ok = False
print('quadrilateral vertices (0,0),(1/2,0),(1,1),(0,2) are exactly the pairwise intersections of the four half-planes')
# source spot check (Section 5.1): path(1, 0.9, 1, 'rlrlrlr') goes negative
t = (Fr(1), Fr(1), Fr(1), Fr(9, 10), Fr(1), Fr(1))
for i in 'rl' + 'rlrlrlr': t = mut(t, i)
print('path(1, 0.9, 1, rlrlrlr) last triple shadows:', [float(t[k]) for k in (1, 3, 5)])
if min(t[k] for k in (1, 3, 5)) >= 0: print('SOURCE_CHECK_FAIL'); ok = False
print('ALL_OK' if ok else 'FAIL')
