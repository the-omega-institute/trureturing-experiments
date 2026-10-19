#!/usr/bin/env python3
"""Samadzadeh-Mojdeh-Nadimi (arXiv:2402.10842v4; AKCE IJGC 22 (2025) 43), Conjecture 1: PC(T(h)) = 5 if h = 2, 3 if h is
odd, 0 otherwise, for the perfect binary tree T(h) (heap labels 1..2^(h+1)-1, children 2v, 2v+1).

A paired dominating set is a dominating set whose induced subgraph has a perfect matching (forest induced subgraphs:
greedy leaf matching decides perfect matchability exactly). For h = 3t - 1 (t = 1, 2, 3) the partition A = union of
{r, 2r, 2r+1} over vertices r at depths 0, 3, ..., 3t-3, B_j = {4r + j} (j = 0..3) is checked to be a pc-partition with
5 parts: no part is paired dominating and every A u B_j is. At h = 5 (odd) this gives PC(T(5)) >= 5 != 3, at h = 8 (even)
PC(T(8)) >= 5 != 0. The search seat's explicit T(5) witness is checked as well. Exit 0 and RESULT PASS when all hold.
"""
import itertools, sys
def tree(h):
    n=2**(h+1)-1; adj={v:set() for v in range(1,n+1)}
    for v in range(1,2**h):
        for c in (2*v,2*v+1): adj[v].add(c); adj[c].add(v)
    return n,adj
def dominating(S,adj,n): return all(v in S or adj[v]&S for v in range(1,n+1))
def has_perfect_matching(S,adj):
    # tree-induced subgraph (forest): greedy leaf matching is exact on forests
    S=set(S); rem=set(S)
    deg={v:len(adj[v]&rem) for v in rem}
    while rem:
        leaves=[v for v in rem if len(adj[v]&rem)<=1]
        if not leaves: return False
        v=leaves[0]
        nb=adj[v]&rem
        if not nb: return False
        u=nb.pop(); rem.discard(v); rem.discard(u)
    return True
def pd(S,adj,n): return dominating(S,adj,n) and has_perfect_matching(S,adj)
def family(t):
    h=3*t-1; n,adj=tree(h)
    depth=lambda v: v.bit_length()-1
    R=[v for v in range(1,n+1) if depth(v)%3==0 and depth(v)<=3*t-3]
    A=set(); B=[set() for _ in range(4)]
    for r in R:
        A|={r,2*r,2*r+1}
        for j in range(4): B[j].add(4*r+j)
    return h,n,adj,[A]+B
ok=True
for t in (1,2,3):
    h,n,adj,parts=family(t)
    allv=set().union(*parts)
    part_ok = allv==set(range(1,n+1)) and sum(len(p) for p in parts)==n and all(parts)
    nonpd = all(not pd(p,adj,n) for p in parts)
    partners = all(any(pd(parts[i]|parts[j],adj,n) for j in range(5) if j!=i) for i in range(5))
    print(f'h={h} n={n} sizes={[len(p) for p in parts]} partition={part_ok} no-part-PD={nonpd} every-part-has-partner={partners}')
    ok = ok and part_ok and nonpd and partners
# explicit T(5) witness from the search seat
n,adj=tree(5)
A={1,2,3,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31}
Bs=[{4,32,36,40,44,48,52,56,60},{5,33,37,41,45,49,53,57,61},{6,34,38,42,46,50,54,58,62},{7,35,39,43,47,51,55,59,63}]
parts=[A]+Bs
w = set().union(*parts)==set(range(1,64)) and sum(map(len,parts))==63 and all(not pd(p,adj,n) for p in parts) and all(pd(A|B,adj,n) for B in Bs)
print('explicit T(5) witness valid:',w); ok=ok and w
print('RESULT PASS' if ok else 'RESULT FAIL'); sys.exit(0 if ok else 1)
