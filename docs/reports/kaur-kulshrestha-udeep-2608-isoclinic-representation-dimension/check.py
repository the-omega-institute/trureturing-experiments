#!/usr/bin/env python3
"""Exact finite-group checks for Kaur--Kulshrestha--Udeep Question (2).

Python >=3.12; SymPy is used only to construct integer cyclotomic polynomials.
No downloaded group tables, GAP output, floating point, or random sampling.
"""
from functools import lru_cache
from itertools import permutations, product
from math import lcm
from sympy import Poly, Symbol, cyclotomic_poly


class Cyclotomic:
    def __init__(self, n):
        self.n = n
        self.poly = tuple(reversed([int(x) for x in Poly(cyclotomic_poly(n, Symbol('x'))).all_coeffs()]))
        self.degree = len(self.poly) - 1
        self.zero = (0,) * self.degree
        self.one = self.integer(1)
        self.roots = tuple(self.reduce([0] * k + [1]) for k in range(n))

    def reduce(self, a):
        a = list(a) + [0] * max(0, self.degree - len(a))
        for i in range(len(a) - 1, self.degree - 1, -1):
            for j in range(self.degree):
                a[i - self.degree + j] -= a[i] * self.poly[j]
        return tuple(a[:self.degree])

    def integer(self, n):
        return (n,) + (0,) * (self.degree - 1)

    def add(self, a, b):
        return tuple(x + y for x, y in zip(a, b))

    @lru_cache(None)
    def mul(self, a, b):
        c = [0] * (2 * self.degree - 1)
        for i, x in enumerate(a):
            for j, y in enumerate(b):
                c[i + j] += x * y
        return self.reduce(c)

    def root(self, k):
        return self.roots[k % self.n]

    def conj(self, a):
        ans = self.zero
        for k, x in enumerate(a):
            ans = self.add(ans, tuple(x * y for y in self.root(-k)))
        return ans

    def embed(self, a, old):
        assert self.n % old.n == 0
        ans = self.zero
        for k, x in enumerate(a):
            ans = self.add(ans, tuple(x * y for y in self.root(k * (self.n // old.n))))
        return ans

    def identity(self, d):
        return tuple(tuple(self.one if i == j else self.zero for j in range(d)) for i in range(d))

    def matmul(self, a, b):
        d = len(a)
        return tuple(tuple(self.sum(self.mul(a[i][k], b[k][j]) for k in range(d))
                           for j in range(d)) for i in range(d))

    def sum(self, xs):
        ans = self.zero
        for x in xs:
            ans = self.add(ans, x)
        return ans


class Group:
    def __init__(self, name, elements, mul, field, reps):
        self.name, self.elements, self.field = name, tuple(elements), field
        idx = {a: i for i, a in enumerate(self.elements)}
        self.n = len(self.elements)
        self.table = tuple(tuple(idx[mul(a, b)] for b in self.elements) for a in self.elements)
        assert self.table[0] == tuple(range(self.n))
        assert all(self.table[i][0] == i for i in range(self.n))
        self.inv = tuple(next(j for j in range(self.n) if self.table[i][j] == 0) for i in range(self.n))
        self.reps = tuple(tuple(r(a) for a in self.elements) for r in reps)
        self.center = tuple(i for i in range(self.n) if all(self.table[i][j] == self.table[j][i]
                                                          for j in range(self.n)))
        self.derived = self.closure(self.comm(i, j) for i in range(self.n) for j in range(self.n))

    def comm(self, i, j):
        t = self.table
        return t[t[t[i][j]][self.inv[i]]][self.inv[j]]

    def closure(self, gens):
        gens = tuple(set(gens))
        seen, queue = {0}, [0]
        for x in queue:
            for g in gens:
                y = self.table[x][g]
                if y not in seen:
                    seen.add(y)
                    queue.append(y)
        return tuple(sorted(seen))

    def validate(self):
        t, f, n = self.table, self.field, self.n
        assert all(t[t[a][b]][c] == t[a][t[b][c]]
                   for a in range(n) for b in range(n) for c in range(n))
        chars, kernels, dims = [], [], []
        for rep in self.reps:
            d = len(rep[0]); ident = f.identity(d)
            assert rep[0] == ident
            assert all(f.matmul(rep[a], rep[b]) == rep[t[a][b]] for a in range(n) for b in range(n))
            char = tuple(f.sum(m[i][i] for i in range(d)) for m in rep)
            chars.append(char)
            kernels.append(sum(1 << a for a in range(n) if rep[a] == ident))
            dims.append(d)
        for i, x in enumerate(chars):
            for j, y in enumerate(chars):
                assert f.sum(f.mul(x[a], f.conj(y[a])) for a in range(n)) == f.integer(n if i == j else 0)
        assert sum(d * d for d in dims) == n
        # Irreducibility, inequivalence and the degree-square sum certify completeness.
        best = {(1 << n) - 1: 0}
        for ker, d in zip(kernels, dims):
            for mask, cost in list(best.items()):
                target = mask & ker
                best[target] = min(best.get(target, n + 1), cost + d)
        assert 1 in best
        self.delta = max(1, best[1])  # source's smallest POSITIVE dimension, including G=1
        self.dims = tuple(dims)
        print(f'GROUP {self.name}: order={n} center={len(self.center)} derived={len(self.derived)} '
              f'irreducible_degrees={dims} delta={self.delta}')


def cyclic(n):
    f = Cyclotomic(n)
    return Group(f'C{n}', range(n), lambda a, b: (a + b) % n, f,
                 [lambda a, k=k: ((f.root(k * a),),) for k in range(n)])


def dihedral(m):
    f = Cyclotomic(lcm(m, 2)); z = f.zero
    es = list(product(range(2), range(m)))
    mul = lambda x, y: ((x[0] + y[0]) % 2, (x[1] + (-1 if x[0] else 1) * y[1]) % m)
    reps = []
    for r in ([0, 1] if m % 2 == 0 else [0]):
        for s in range(2):
            reps.append(lambda x, r=r, s=s: ((f.root((f.n // 2) * (s*x[0]+r*x[1])),),))
    for k in range(1, (m - 1)//2 + 1):
        def rho(x, k=k):
            a, b = f.root((f.n // m) * k * x[1]), f.root(-(f.n // m) * k * x[1])
            return ((a, z), (z, b)) if x[0] == 0 else ((z, a), (b, z))
        reps.append(rho)
    return Group(f'D{2*m}', es, mul, f, reps)


def dicyclic(m):
    f = Cyclotomic(lcm(2*m, 4)); z = f.zero
    es = list(product(range(2), range(2*m)))
    mul = lambda x, y: ((x[0] + y[0]) % 2,
                        (x[1]+(-1 if x[0] else 1)*y[1]+m*x[0]*y[0]) % (2*m))
    reps = []
    for r in range(2):
        for s in range(4):
            if (2*s - 2*r*m) % 4 == 0:
                reps.append(lambda x, r=r, s=s: ((f.root((f.n//2)*r*x[1]+(f.n//4)*s*x[0]),),))
    for k in range(1, m):
        def rho(x, k=k):
            a, b = f.root((f.n//(2*m))*k*x[1]), f.root(-(f.n//(2*m))*k*x[1])
            return ((a, z), (z, b)) if x[0] == 0 else ((z, f.mul(f.integer((-1)**k), a)), (b, z))
        reps.append(rho)
    return Group('Q8' if m == 2 else f'Dic{4*m}', es, mul, f, reps)


def direct(a, b):
    f = Cyclotomic(lcm(a.field.n, b.field.n))
    def tensor(x, y):
        return tuple(tuple(f.mul(f.embed(x[i][j], a.field), f.embed(y[k][l], b.field))
                           for j in range(len(x)) for l in range(len(y)))
                     for i in range(len(x)) for k in range(len(y)))
    reps = [lambda x, r=r, s=s: tensor(r[x[0]], s[x[1]]) for r in a.reps for s in b.reps]
    return Group(a.name+'x'+b.name, product(range(a.n), range(b.n)),
                 lambda x, y: (a.table[x[0]][y[0]], b.table[x[1]][y[1]]), f, reps)


def alternating_four():
    es = [p for p in permutations(range(4)) if sum(p[i] > p[j] for i in range(4) for j in range(i+1,4)) % 2 == 0]
    compose = lambda a, b: tuple(a[b[i]] for i in range(4))
    f = Cyclotomic(3)
    klein = [p for p in es if compose(p, p) == es[0]]
    c = (1, 2, 0, 3)
    cosets = [{compose(v, r) for v in klein} for r in [es[0], c, compose(c, c)]]
    power = {x: k for k, coset in enumerate(cosets) for x in coset}
    reps = [lambda x, k=k: ((f.root(k*power[x]),),) for k in range(3)]
    # Action on the basis e_i-e_3 (i=0,1,2) of the sum-zero permutation module.
    reps.append(lambda x: tuple(tuple(f.integer(int(x[j] == i)-int(x[3] == i))
                                    for j in range(3)) for i in range(3)))
    return Group('A4', es, compose, f, reps)


def subgroup_table(g, subgroup):
    idx = {x: i for i, x in enumerate(subgroup)}
    return tuple(tuple(idx[g.table[x][y]] for y in subgroup) for x in subgroup)


def quotient(g):
    cosets, at = [], {}
    for x in range(g.n):
        if x not in at:
            c = {g.table[x][z] for z in g.center}
            for y in c:
                at[y] = len(cosets)
            cosets.append(min(c))
    table = tuple(tuple(at[g.table[x][y]] for y in cosets) for x in cosets)
    return table, cosets


def isomorphisms(a, b):
    if len(a) != len(b):
        return
    n = len(a)
    def order(t, x):
        y = x; k = 1
        while y:
            y = t[y][x]; k += 1
        return k
    oa, ob = [order(a, i) for i in range(n)], [order(b, i) for i in range(n)]
    def grow(mapping):
        if len(mapping) == n:
            yield tuple(mapping[i] for i in range(n)); return
        x = next(i for i in range(n) if i not in mapping)
        for y in range(n):
            if y in mapping.values() or oa[x] != ob[y]:
                continue
            trial = dict(mapping); trial[x] = y
            changed, good = True, True
            while changed and good:
                changed = False
                for u, v in product(tuple(trial), repeat=2):
                    uv, image = a[u][v], b[trial[u]][trial[v]]
                    if uv in trial:
                        if trial[uv] != image:
                            good = False; break
                    elif image in trial.values() or oa[uv] != ob[image]:
                        good = False; break
                    else:
                        trial[uv] = image; changed = True
            if good:
                yield from grow(trial)
    yield from grow({0: 0})


def isoclinic(g, h):
    if len(g.derived) != len(h.derived):
        return False
    qg, rg = quotient(g); qh, rh = quotient(h)
    for phi in isomorphisms(qg, qh):
        generators = {(g.comm(x, y), h.comm(rh[phi[i]], rh[phi[j]]))
                      for i, x in enumerate(rg) for j, y in enumerate(rg)}
        seen, queue = {(0, 0)}, [(0, 0)]
        for a, b in queue:
            for c, d in generators:
                v = (g.table[a][c], h.table[b][d])
                if v not in seen:
                    seen.add(v); queue.append(v)
        if len(seen) == len(g.derived) == len({x for x,y in seen}) == len({y for x,y in seen}):
            return True
    return False


def main():
    cyclics = {n: cyclic(n) for n in range(1, 16)}
    c2, c3, c4 = [cyclics[n] for n in [2, 3, 4]]
    d6, d8, q8, d12, dic12 = dihedral(3), dihedral(4), dicyclic(2), dihedral(6), dicyclic(3)
    # Standard complete isomorphism list of the 28 groups of orders 1 through 15.
    small = list(cyclics.values()) + [direct(c2,c2), direct(c4,c2), direct(direct(c2,c2),c2),
                                    direct(c3,c3), direct(cyclics[6],c2),
                                    d6,d8,q8,dihedral(5),alternating_four(),d12,dic12,dihedral(7)]
    assert len(small) == 28
    extra = [direct(g,a) for g,a in [(d8,c2),(q8,c2),(d8,c3),(q8,c3),
                                   (d12,c2),(dic12,c2),(d8,c4),(q8,c4)]]
    fidelity = direct(d6,c2)
    groups = small + extra + [fidelity]
    for g in groups:
        g.validate()
    # Section 3, Example on failure of product additivity: delta(S3 x C2)=2,
    # while delta(S3)+delta(C2)=3. Also the Section 1 C2 x C2 example has no
    # faithful irreducible representation.
    assert fidelity.delta == 2 and d6.delta + c2.delta == 3
    v4 = next(g for g in small if g.name == 'C2xC2')
    assert v4.delta == 2 and all(len(r[0]) == 1 and any(r[x] == r[0] for x in range(1,4)) for r in v4.reps)
    print('FIDELITY: delta(S3 x C2)=2; delta(S3)+delta(C2)=3; C2xC2 has no faithful irreducible')
    eligible = 0; nonidentical = []
    for i, g in enumerate(groups):
        for h in groups[i:]:
            if g.n != h.n:
                continue
            if next(isomorphisms(subgroup_table(g,g.center), subgroup_table(h,h.center)), None) is None:
                continue
            if isoclinic(g, h):
                assert g.delta == h.delta, (g.name,h.name,g.delta,h.delta)
                eligible += 1
                if g is not h:
                    nonidentical.append((g.name,h.name,g.delta))
    print(f'IDENTITY: all {eligible} eligible unordered pairs (self-pairs included) in {len(groups)} presentations')
    for a,b,d in nonidentical:
        print(f'PAIR {a} ~ {b}: equal centers, delta={d}')
    assert ('D8','Q8',2) in nonidentical and ('D12','Dic12',2) in nonidentical
    print('COVERAGE: all 28 isomorphism types with order <=15; eight central direct-factor cases and S3xC2')
    print('LIMITATION: finite sanity checks do not prove the universal claim')
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
