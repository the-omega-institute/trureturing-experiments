#!/usr/bin/env python3
"""Exact finite checks for Clarke, arXiv:2610.17988v1, Conjecture 8.5.

No numerical hull library or floating-point arithmetic is used. Supporting
facets are enumerated using SymPy rational nullspaces; integer inequalities
then decide lattice membership. Finite checks do not prove the conjecture.
"""
from itertools import combinations, product
from math import comb, gcd, lcm
from functools import reduce
from random import Random
import sympy as sp

SEED = 261085


def convolve(a, b, K):
    return [sum(a[j] * b[i-j] for j in range(i+1)
                if j < len(a) and i-j < len(b)) for i in range(K+1)]


def factor(sign, power, stride=1):
    a = [0] * (stride * power + 1)
    for j in range(power+1):
        a[stride*j] = comb(power, j) * sign**j
    return a


def radd(a, b):
    return (a[0]+b[0], a[1]+b[1])


def rmul(a, b):
    return (a[0]*b[0]+a[1]*b[1], a[0]*b[1]+a[1]*b[0])


def rconv(a, b, K):
    out = []
    for i in range(K+1):
        v = (0, 0)
        for j in range(i+1):
            if j < len(a) and i-j < len(b):
                v = radd(v, rmul(a[j], b[i-j]))
        out.append(v)
    return out


class Polytope:
    def __init__(self, vertices, A, b):
        self.V = sorted(set(tuple(v) for v in vertices))
        self.d = len(self.V[0])
        self.A, self.b = sp.Matrix(A), tuple(b)
        assert self.A * self.A == sp.eye(self.d)
        assert self.A * sp.Matrix(self.b) + sp.Matrix(self.b) == sp.zeros(self.d, 1)
        assert self.A != sp.eye(self.d) or any(self.b)
        assert sp.Matrix([[v[i]-self.V[0][i] for i in range(self.d)]
                          for v in self.V[1:]]).rank() == self.d
        assert {self.act(v, 1) for v in self.V} == set(self.V)
        self.facets = self.get_facets()
        self.bounds = [(min(v[i] for v in self.V), max(v[i] for v in self.V))
                       for i in range(self.d)]
        homogeneous = self.A.row_join(sp.Matrix(self.b)).col_join(
            sp.Matrix([[0]*self.d + [1]]))
        self.q = int((homogeneous-sp.eye(self.d+1)).rank())
        self.p = self.d + 1 - self.q
        assert self.p >= 1 and self.q >= 1
        self.cache = {}

    def act(self, x, k):
        return tuple(sum(int(self.A[i,j])*x[j] for j in range(self.d))
                     + k*self.b[i] for i in range(self.d))

    def get_facets(self):
        facets = set()
        for S in combinations(self.V, self.d):
            M = sp.Matrix([list(v)+[1] for v in S])
            null = M.nullspace()
            if len(null) != 1:
                continue
            v = null[0]
            den = lcm(*(int(x.q) for x in v))
            a = [int(x*den) for x in v]
            g = reduce(gcd, (abs(x) for x in a))
            a = [x//g for x in a]
            values = [sum(x*y for x,y in zip(a, list(w)+[1])) for w in self.V]
            if max(values) <= 0:
                facets.add(tuple(a))
            elif min(values) >= 0:
                facets.add(tuple(-x for x in a))
        assert len(facets) >= self.d + 1
        return sorted(facets)

    def points(self, k):
        if k not in self.cache:
            if k == 0:
                self.cache[k] = {(0,)*self.d}
            else:
                box = product(*(range(k*lo, k*hi+1) for lo,hi in self.bounds))
                self.cache[k] = {x for x in box if all(
                    sum(a[i]*x[i] for i in range(self.d)) + k*a[-1] <= 0
                    for a in self.facets)}
        return self.cache[k]

    def counts(self, K):
        L, N = [], []
        for k in range(K+1):
            pts = self.points(k)
            assert {self.act(x,k) for x in pts} == pts
            L.append(len(pts))
            N.append(sum(self.act(x,k) == x for x in pts))
            assert (L[-1]-N[-1]) % 2 == 0
        return L, N

    def phi(self, K):
        L, N = self.counts(K)
        E = [((a+b)//2, (a-b)//2) for a,b in zip(L,N)]
        D = [(1,0)]
        for _ in range(self.p):
            D = rconv(D, [(1,0),(-1,0)], len(D))
        for _ in range(self.q):
            D = rconv(D, [(1,0),(0,-1)], len(D))
        return rconv(D,E,K)

    def numerators(self, K):
        L,N = self.counts(K)
        return (convolve(factor(-1,self.d+1),L,K),
                convolve(factor(-1,self.p,2),N,K))


def check_bound(P, K):
    h,H = P.numerators(K)
    assert min(h) >= 0 and min(H) >= 0
    if P.p >= P.q:
        left,right = convolve(factor(1,P.p-P.q),h,K), H
    else:
        left,right = h, convolve(factor(1,P.q-P.p),H,K)
    assert all(a >= b for a,b in zip(left,right)), (P.V,left,right)
    phi = P.phi(K)
    assert [a+b for a,b in phi] == h
    assert convolve(factor(1,P.p),[a-b for a,b in phi],K) == convolve(factor(1,P.q),H,K)
    return h,H


def bipyramid(P, v):
    w = P.act(v,1)
    vertices = [x+(0,) for x in P.V] + [v+(1,), w+(-1,)]
    A = sp.diag(P.A, sp.Matrix([[-1]]))
    return Polytope(vertices,A,list(P.b)+[0])


def check_construction(P, K=4):
    v = P.V[0]
    B = bipyramid(P,v)
    assert (B.p,B.q) == (P.p,P.q+1)
    for k in range(K+1):
        points = B.points(k)
        assert {z for z in points if z[-1] == 0} == {x+(0,) for x in P.points(k)}
        for j in range(1,k+1):
            for sign,apex in [(1,v),(-1,P.act(v,1))]:
                expected = {tuple(x[i]+j*apex[i] for i in range(P.d))+(sign*j,)
                            for x in P.points(k-j)}
                assert {z for z in points if z[-1] == sign*j} == expected
    LP,NP = P.counts(K)
    LB,NB = B.counts(K)
    assert NP == NB
    assert LB == [LP[k]+2*sum(LP[:k]) for k in range(K+1)]
    assert B.phi(K) == rconv([(1,0),(1,0)],P.phi(K),K)
    return B


def trimmed(a):
    a = list(a)
    while len(a)>1 and a[-1] == 0:
        a.pop()
    return a


def source_checks():
    square = list(product(range(2),repeat=2))
    orthogonal = Polytope(square,[[-1,0],[0,1]],[1,0])
    diagonal = Polytope(square,[[0,1],[1,0]],[0,0])
    pyramid = Polytope([(0,0,0),(1,0,0),(0,1,0),(1,0,1),(0,1,1)],
                       [[0,1,0],[1,0,0],[0,0,1]],[0,0,0])
    source_bipyramid = Polytope([x+(0,) for x in square]+[(0,0,1),(1,1,-1)],
                               [[-1,0,-1],[0,1,0],[0,0,1]],[1,0,0])
    sources = [('orthogonal square',orthogonal,[1,1],[1,0,1]),
               ('diagonal square',diagonal,[1,1],[1,2,1]),
               ('pyramid',pyramid,[1,1],[1,1,1,1]),
               ('bipyramid',source_bipyramid,[1,2,1],[1,2,2,2,1])]
    for name,P,h_expected,H_expected in sources:
        K = 10
        h,H = check_bound(P,K)
        assert trimmed(h) == h_expected, (name,h)
        assert trimmed(H) == H_expected, (name,H)
        if name == 'pyramid':
            assert P.counts(K)[1] == [(k//2+1)**2 for k in range(K+1)]
        print('SOURCE',name,'p,q=',P.p,P.q,'L,N=',P.counts(4),
              'h*=',trimmed(h),'H=',trimmed(H),'through degree',K)
    # Example 5.6 and its displayed Conjecture 8.5 comparison.
    W = [(1,0,0,0,0,0),(0,0,1,0,0,0),
         (0,1,0,1,1,0),(0,0,1,1,1,0)]
    swap = lambda x:x[3:]+x[:3]
    P = Polytope(W+[swap(x) for x in W],
                 [[int(j == (i+3)%6) for j in range(6)] for i in range(6)],[0]*6)
    h,H = check_bound(P,6)
    assert trimmed(h) == [1,1,0,1,1]
    assert trimmed(H) == [1,0,1,1,0,1]
    assert convolve([1,1],h,6) == [1,2,1,1,2,1,0]
    assert P.phi(6)[2] == (1,-1)
    print('SOURCE Example 5.6 p,q=',P.p,P.q,'L,N=',P.counts(6),
          'h*=',trimmed(h),'H=',trimmed(H),'phi_2=',P.phi(6)[2])
    return [x[1] for x in sources]


def small_checks(sources):
    # Every reflected integer interval with endpoints in [-3,3].
    intervals = []
    for a in range(-3,4):
        for b in range(a+1,4):
            intervals.append(Polytope([(a,),(b,)],[[-1]],[a+b]))
    # Every full-dimensional invariant subset of the 3x3 integer grid,
    # for reflection in x=0 and diagonal interchange. Duplicated hulls
    # remain independent generating-set checks.
    grid = list(product(range(-1,2),repeat=2))
    planar = []
    for A in [sp.diag(-1,1),sp.Matrix([[0,1],[1,0]])]:
        unseen = set(grid)
        orbits = []
        while unseen:
            x = min(unseen)
            y = tuple(int(z) for z in A*sp.Matrix(x))
            orb = {x,y}
            unseen -= orb
            orbits.append(orb)
        for mask in range(1,1<<len(orbits)):
            V = sorted(set().union(*(orbits[i] for i in range(len(orbits)) if mask>>i&1)))
            if len(V) < 3 or sp.Matrix([[x[i]-V[0][i] for i in range(2)] for x in V[1:]]).rank()!=2:
                continue
            planar.append(Polytope(V,A,[0,0]))
    # 24 deterministic random full-dimensional 3-polytopes, including
    # affine reflections with no fixed lattice point at height one.
    rng = Random(SEED)
    spatial = []
    actions = [(sp.diag(-1,1,1),[1,0,0]),
               (sp.diag(-1,-1,1),[0,0,0]),
               (-sp.eye(3),[0,0,0]),
               (sp.Matrix([[0,1,0],[1,0,0],[0,0,1]]),[0,0,0])]
    while len(spatial)<24:
        A,b = actions[len(spatial)%len(actions)]
        seeds = [tuple(rng.randrange(-1,2) for _ in range(3)) for _ in range(4)]
        V = sorted(set(seeds+[tuple(int(z) for z in A*sp.Matrix(x)+sp.Matrix(b)) for x in seeds]))
        if sp.Matrix([[x[i]-V[0][i] for i in range(3)] for x in V[1:]]).rank()!=3:
            continue
        spatial.append(Polytope(V,A,b))
    all_cases = intervals+planar+spatial
    branches = {'p>=q':0,'p<q':0}
    for P in all_cases:
        check_bound(P,2*(P.d+1)+2)
        branches['p>=q' if P.p>=P.q else 'p<q'] += 1
    selected = sources+intervals[:3]+planar[::max(1,len(planar)//10)]+spatial
    constructions = 0
    for P in selected:
        B = check_construction(P)
        constructions += 1
        # Full balancing iteration is directly re-enumerated through k=2.
        for _ in range(1,max(P.p-P.q,0)):
            B = check_construction(B,K=2)
            constructions += 1
        assert B.p <= B.q or P.p <= P.q
        assert all(a>=0 and b>=0 for a,b in B.phi(2))
    print('SMALL exhaustive intervals=',len(intervals),'planar generating sets=',len(planar),
          'random spatial=',len(spatial),'seed=',SEED,'branches=',branches)
    print('CONSTRUCTION directly enumerated steps=',constructions,
          'first steps k=0..4; additional balancing steps k=0..2')
    print('All small input coefficient windows i=0..2(d+1)+2 passed; finite evidence only.')


def main():
    # Negative controls: ordinary effectiveness does fail; coefficient
    # domination cannot be cancelled through a positive factor.
    assert convolve([1,1],[1,-1,1],3) == [1,0,0,1]
    assert min([1,-1,1]) < 0
    sources = source_checks()
    small_checks(sources)
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
