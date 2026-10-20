"""Shala, arXiv:2112.13766v6, Question 7.4: an atomistic lattice that is weakly but not strongly coset-like.

For m >= 4, L_m = {empty, X, A, Omega} together with the 2m singletons of Omega = {1..2m}, with X = {1..m-1} and
A = {1..m}, ordered by inclusion. Checks for m = 4..8: L_m is a lattice (closed under intersection, unique least upper
bounds), atomistic, its join-irreducibles are exactly the atoms, the local sums S(q) of P(L_m, s) vanish for every
non-integer q (so P(L_m, s) = 1 - 2^{-s} - m (2m)^{-s} is a finite ordinary Dirichlet series), and |J_X| = m - 1 does not
divide |J| = 2m.
"""
from fractions import Fraction
import itertools
def build(m):
    O=frozenset(range(1,2*m+1)); X=frozenset(range(1,m)); A=frozenset(range(1,m+1))
    L=[frozenset(),X,A,O]+[frozenset([i]) for i in range(1,2*m+1)]
    return L,O
def check(m):
    L,O=build(m); S=set(L)
    # lattice: closed under intersection and has top; join = smallest member containing union
    closed=all((a&b) in S for a in L for b in L)
    def join(a,b): return min([c for c in L if a|b<=c], key=len)
    # uniqueness of least upper bound
    lub_ok=all(len([c for c in L if a|b<=c and len(c)==len(join(a,b))])==1 for a in L for b in L)
    bot=frozenset()
    atoms=[x for x in L if x!=bot and all(not (bot<y<x) for y in L)]
    def supirred(x): return x!=bot and all(not(join(a,b)==x and a!=x and b!=x) for a in L for b in L)
    J=[x for x in L if supirred(x)]
    atomistic=all(x==bot or x==frozenset().union(*[a for a in atoms if a<=x]) and join_all([a for a in atoms if a<=x],join,bot)==x for x in L)
    # Möbius mu(x, top)
    mu={}
    for x in sorted(L,key=lambda t:-len(t)):
        mu[x]=1 if x==O else -sum(mu[y] for y in L if x<y and y<=O)
    Jx=lambda x:[j for j in J if j<=x]
    S_q={}
    for x in L:
        if x==bot: continue
        q=Fraction(len(J),len(Jx(x)))
        S_q[q]=S_q.get(q,0)+mu[x]
    weak=all(v==0 for q,v in S_q.items() if q.denominator!=1)
    strong=all(len(J)%len(Jx(x))==0 for x in L if x!=bot)
    return closed,lub_ok,set(J)==set(atoms),atomistic,{str(k):v for k,v in S_q.items() if v!=0},weak,strong,len(L)
def join_all(xs,join,bot):
    r=bot
    for x in xs: r=join(r,x)
    return r

import sys
ok = True
for m in range(4, 9):
    closed, lub_ok, j_atoms, atomistic, sq, weak, strong, size = check(m)
    expected = {'1': 1, '2': -1, str(2 * m): -m}
    good = closed and lub_ok and j_atoms and atomistic and weak and not strong and sq == expected
    ok = ok and good
    print(f"m={m}: |L|={size}, lattice {closed and lub_ok}, atomistic {atomistic}, J = atoms {j_atoms}, nonzero S(q) {sq}, weakly {weak}, strongly {strong}")
print("RESULT PASS" if ok else "RESULT FAIL")
sys.exit(0 if ok else 1)
