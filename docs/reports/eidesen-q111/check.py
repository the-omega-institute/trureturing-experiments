#!/usr/bin/env python3
"""Exact monomial checks for Eidesen Q11.1; Python standard library only.

M=(p,e) means M e_j = z_N**e[j] e_{p[j]}, N=4m, z_N primitive.
All equality tests are integer congruences. No floating-point arithmetic is used.
A single nonzero entry per column makes unitarity automatic after permutation
and exponent validation. Distinct diagonal weights and connected permutation
graphs certify scalar commutants, hence irreducibility for these unitary blocks.
The full model at m=5 is enumerated independently of the product count formula.
"""
from itertools import product
from math import gcd
import json


def identity(d):
    return tuple(range(d)), (0,) * d


def mul(a, b, n):
    p, e = a
    q, f = b
    return tuple(p[q[j]] for j in range(len(q))), tuple((f[j] + e[q[j]]) % n for j in range(len(q)))


def power(a, k, n):
    v = identity(len(a[0]))
    for _ in range(k):
        v = mul(v, a, n)
    return v


def projkey(a, n):
    p, e = a
    return p, tuple((x - e[0]) % n for x in e)


def tensor(a, b, n):
    p, e = a
    q, f = b
    return tuple(i * len(q) + j for i in p for j in q), tuple((i + j) % n for i in e for j in f)


def islogical(a, w):
    return {a[0][j] for j in w} == set(w)


def scalar_phase(a, w):
    assert w
    p, e = a
    if all(p[j] == j and e[j] == e[w[0]] for j in w):
        return e[w[0]]
    return None


def commutant_one(diagonal, generators, support):
    """Distinct diagonal eigenvalues force diagonal commutants; edges equate entries."""
    p, e = diagonal
    assert all(p[j] == j for j in support)
    assert len({e[j] for j in support}) == len(support)
    seen = {support[0]}
    while True:
        new = seen | {g[0][j] for g in generators for j in seen}
        if new == seen:
            break
        seen = new
    assert seen == set(support)


def invariants(matrices, w):
    logical = [i for i, a in enumerate(matrices) if islogical(a, w)]
    scalar = [(i, scalar_phase(a, w)) for i, a in enumerate(matrices) if scalar_phase(a, w) is not None]
    d = len(matrices[0][0])
    eigencoords = [j for j in range(d) if all(matrices[i][0][j] == j and matrices[i][1][j] == c for i, c in scalar)]
    # These scalar matrices must be diagonal on the WHOLE ambient space before
    # counting common eigenspace dimension by surviving coordinate vectors.
    assert all(matrices[i][0] == tuple(range(d)) for i, _ in scalar)
    return {'L': len(logical), 'S': len(scalar), 'E': len(eigencoords)}


def clifford_partition(matrices, w, diagonal, generators):
    stabilizer = [a for a in matrices if islogical(a, w)]
    assert stabilizer
    commutant_one(diagonal, generators, w)
    images = {tuple(sorted(a[0][j] for j in w)) for a in matrices}
    assert len(images) * len(stabilizer) == len(matrices)
    assert sum(map(len, images)) == len(matrices[0][0])
    assert set().union(*map(set, images)) == set(range(len(matrices[0][0])))
    # Orthogonal coset images identify the literal induced representation with
    # the ambient one, via delta_g tensor w |-> pi(g)w, for these explicit blocks.


def build(m, r):
    n = 4 * m
    X = ((1, 0), (0, 0))
    P = ((0, 1), (0, m))
    C = ((2, 3, 0, 1), (0, 0, 0, 0))
    XX = ((1, 0, 3, 2), (0, 0, 0, 0))
    PP = ((0, 1, 2, 3), (0, m, 2*m, 3*m))
    Bdiag = ((0, 1), (4, (4*(m-1)) % n))
    Fdiag = (tuple(range(4)), tuple((4*pow(r, j, m)) % n for j in range(4)))
    Fshift = ((3, 0, 1, 2), (0, 0, 0, 0))
    # Elements in normal form c^k b^l a^u for A, b^k a^u otherwise.
    Aelts = list(product(range(2), range(2), range(4)))
    Qelts = list(product(range(2), range(4)))
    Belts = list(product(range(2), range(m)))
    Felts = list(product(range(4), range(m)))
    A = [mul(mul(power(C, k, n), power(XX, l, n), n), power(PP, u, n), n) for k,l,u in Aelts]
    Q = [mul(power(X, k, n), power(P, u, n), n) for k,u in Qelts]
    B = [mul(power(X, k, n), power(Bdiag, u, n), n) for k,u in Belts]
    F = [mul(power(Fshift, k, n), power(Fdiag, u, n), n) for k,u in Felts]
    def dihmul(x, y, order):
        k,u = x; l,v = y
        return (k+l) % 2, (((-1)**l)*u+v) % order
    def amul(x,y):
        k,l,u=x; h,i,v=y
        return ((k+h)%2, *dihmul((l,u),(i,v),4))
    def fmul(x,y):
        k,u=x; l,v=y
        return (k+l)%4, (pow(r, (-l)%4, m)*u+v)%m
    specs = [
        ('A', Aelts,A,amul,PP,[C,XX],(0,1),[XX]),
        ('B', Belts,B,lambda x,y:dihmul(x,y,m),Bdiag,[X],(0,),[]),
        ('Q', Qelts,Q,lambda x,y:dihmul(x,y,4),P,[X],(0,1),[X]),
        ('F', Felts,F,fmul,Fdiag,[Fshift],(0,),[])]
    for name,elts,mats,gmul,diag,gens,w,restricted_gens in specs:
        assert len(mats) == len(set(elts))
        assert len({projkey(a,n) for a in mats}) == len(mats)
        byelt = dict(zip(elts,mats))
        for a in mats:
            assert sorted(a[0]) == list(range(len(a[0])))
            assert all(0 <= e < n for e in a[1])
        for x, a in byelt.items():
            for y, b in byelt.items():
                assert projkey(mul(a,b,n),n) == projkey(byelt[gmul(x,y)],n), (name,x,y)
        commutant_one(diag,gens,tuple(range(len(diag[0]))))
        clifford_partition(mats,w,diag,restricted_gens)
    assert all(mul(mul(Fshift,Fdiag,n),power(Fshift,3,n),n)[1][j] == power(Fdiag,r,n)[1][j] for j in range(4))
    return n,A,B,Q,F


def full_check(n,A,B,Q,F):
    good = tuple(a*16+b*8+q*4 for a in range(4) for b in range(2) for q in range(2))
    bad = tuple(a*16+q*4+f for a in range(2) for q in range(2) for f in range(4))
    assert len(good)==len(bad)==16
    counts = [[0,0],[0,0]]
    scalar_mats=[[],[]]
    full_keys=set()
    ab=[tensor(a,b,n) for a,b in product(A,B)]
    qf=[tensor(q,f,n) for q,f in product(Q,F)]
    for x,y in product(ab,qf):
        z=tensor(x,y,n)
        full_keys.add(projkey(z,n))
        for i,w in enumerate((good,bad)):
            counts[i][0]+=islogical(z,w)
            c=scalar_phase(z,w)
            if c is not None:
                counts[i][1]+=1
                assert z[0] == tuple(range(64))
                scalar_mats[i].append((z,c))
    assert len(full_keys)==len(A)*len(B)*len(Q)*len(F)
    dims=[sum(all(z[1][j]==c for z,c in ss) for j in range(64)) for ss in scalar_mats]
    return {'G':len(full_keys),'good':{'L':counts[0][0],'S':counts[0][1],'E':dims[0]},'bad':{'L':counts[1][0],'S':counts[1][1],'E':dims[1]}}


def run():
    for t in (1,2,3,4):
        r=2*t; m=r*r+1
        assert m%2==1 and r*r%m==m-1 and pow(r,4,m)==1
        assert len({pow(r,j,m) for j in range(4)})==4 and gcd(r-1,m)==1
        n,A,B,Q,F=build(m,r)
        aa=invariants(A,(0,1)); bb=invariants(B,(0,)); qq=invariants(Q,(0,1)); ff=invariants(F,(0,))
        assert aa=={'L':8,'S':1,'E':4}
        assert bb=={'L':m,'S':m,'E':1}
        assert qq=={'L':8,'S':1,'E':2}
        assert ff=={'L':m,'S':m,'E':1}
        af=invariants(A,(0,1,2,3)); bf=invariants(B,(0,1)); ffull=invariants(F,(0,1,2,3))
        good={'L':af['L']*bf['L']*qq['L']*ff['L'],'S':af['S']*bf['S']*qq['S']*ff['S'],'E':4*2*2*ff['E']}
        bad={'L':aa['L']*bb['L']*qq['L']*ffull['L'],'S':aa['S']*bb['S']*qq['S']*ffull['S'],'E':aa['E']*bb['E']*2*4}
        assert good=={'L':256*m*m,'S':m,'E':16}
        assert bad=={'L':256*m*m,'S':m,'E':32}
        row={'t':t,'m':m,'G':len(A)*len(B)*len(Q)*len(F),'good':good,'bad':bad,'blocks_exact_projective_and_irreducible':True,'clifford_coset_partition_checked':True}
        if t==1:
            direct=full_check(n,A,B,Q,F)
            assert direct=={'G':1024*m*m,'good':good,'bad':bad}
            row['direct_full_64d_enumeration']=direct
            print('FIDELITY Proposition 8.1 at n=2: A has order 16, dim 4; U_A has dim 2, L=8, S=1, maximal eigenspace dim 4',flush=True)
        print(json.dumps(row,sort_keys=True),flush=True)
    print('ALL_CHECKS_OK')


if __name__=='__main__':
    run()
