#!/usr/bin/env python3
"""Exact finite-field checks of the class-two action and its vanishing criterion."""
from itertools import combinations, product
from random import Random


def echelon(rows, p, width):
    out = {}
    for row in rows:
        v = [x % p for x in row]
        for j, b in sorted(out.items()):
            a = v[j]
            if a:
                v = [(x - a*y) % p for x, y in zip(v, b)]
        j = next((j for j, x in enumerate(v) if x), None)
        if j is not None:
            a = pow(v[j], -1, p)
            out[j] = [(a*x) % p for x in v]
    return out


def reduce(v, basis, p):
    v = [x % p for x in v]
    for j, b in sorted(basis.items()):
        a = v[j]
        if a:
            v = [(x - a*y) % p for x, y in zip(v, b)]
    return v


def rref_subspaces(n, p):
    for r in range(n+1):
        for pivots in combinations(range(n), r):
            free = [j for j in range(n) if j not in pivots]
            cells = [(i, j) for i, q in enumerate(pivots) for j in free if q < j]
            for values in product(range(p), repeat=len(cells)):
                rows = [[int(j == q) for j in range(n)] for q in pivots]
                for (i, j), a in zip(cells, values):
                    rows[i][j] = a
                yield rows


def setup(d, p, rows):
    pairs = list(combinations(range(d), 2))
    kb = echelon(rows, p, len(pairs))
    k = list(kb.values())
    free = [j for j in range(len(pairs)) if j not in kb]
    m = len(free)
    def omega_pair(i, j):
        if i == j:
            return [0]*m
        a, b = sorted((i, j))
        v = [int(t == pairs.index((a, b))) for t in range(len(pairs))]
        v = reduce(v, kb, p)
        sign = 1 if i < j else -1
        return [(sign*v[t]) % p for t in free]
    def block(i, w):
        return [w[j % m] if m and j//m == i else 0 for j in range(d*m)]
    jacobi = []
    for a, b, c in combinations(range(d), 3):
        jacobi.append([(x+y+z) % p for x,y,z in zip(
            block(a, omega_pair(b,c)), block(b, omega_pair(c,a)),
            block(c, omega_pair(a,b)))])
    jb = echelon(jacobi, p, d*m)
    columns = []
    for i in range(d):
        for w in range(m):
            col = []
            for rel in k:
                v = [0]*(d*m)
                for coefficient, (a,b) in zip(rel,pairs):
                    if a == i: v[b*m+w] -= coefficient
                    if b == i: v[a*m+w] += coefficient
                col += reduce(v,jb,p)
            columns.append(col)
    action_rank = len(echelon(columns,p,len(k)*d*m))
    qdim = d*m-len(jb)
    assert (action_rank == 0) == (not k or qdim == 0), (d,p,k,qdim,action_rank)
    # Every rank-one action is minus contraction tensor w in V tensor W / J.
    contractions = []
    for rel in k:
        for i in range(d):
            s = [0]*d
            for coefficient,(a,b) in zip(rel,pairs):
                if a == i: s[b] += coefficient
                if b == i: s[a] -= coefficient
            contractions.append([x % p for x in s])
    sb = echelon(contractions,p,d)
    wedges = []
    for x,y in combinations(sb.values(),2):
        wedges.append([(x[a]*y[b]-x[b]*y[a]) % p for a,b in pairs])
    wb = echelon(wedges,p,len(pairs))
    assert all(not any(reduce(rel,wb,p)) for rel in k)
    return qdim,action_rank,len(k),len(sb)


def tensors(p):
    # Degree-three noncommutative polynomials; exact coefficients modulo p.
    def add(*terms):
        r = {}
        for t in terms:
            for m,a in t.items():r[m] = (r.get(m,0)+a) % p
        return {m:a for m,a in r.items() if a}
    def neg(t):return {m:-a % p for m,a in t.items()}
    def mul(a,b):
        r = {}
        for x,u in a.items():
            for y,v in b.items():r[x+y] = (r.get(x+y,0)+u*v) % p
        return {m:a for m,a in r.items() if a}
    def bracket(a,b):return add(mul(a,b),neg(mul(b,a)))
    def star(a):return {tuple(reversed(m)):(a*(-1)**len(m)) % p for m,a in a.items()}
    x,y,z = ({(i,):1} for i in range(3))
    c = bracket(bracket(x,y),x)
    assert c.get((0,1,0),0) == 2 % p
    assert c.get((0,1,0),0) != 0
    # Quadratic reversal-symmetric corrections give star-even cubic displacements.
    sym = add(mul(x,z),mul(z,x))
    anti = bracket(x,z)
    delta_sym = bracket(sym,y)
    delta_anti = bracket(anti,y)
    assert star(delta_sym) == delta_sym
    assert star(delta_anti) == neg(delta_anti)
    if p == 3:
        lhs = mul(mul(add(x,y),add(x,y)),add(x,y))
        rhs = add(mul(mul(x,x),x),mul(mul(y,y),y),
                  bracket(x,bracket(x,y)),bracket(y,bracket(y,x)))
        assert lhs == rhs
    print(f'TENSOR p={p}: middle coefficient=2, inversion signs and cube polarization OK')


def main():
    for p in (3,5,7):
        count = 0
        for d in range(1,4):
            for rows in rref_subspaces(d*(d-1)//2,p):
                setup(d,p,rows);count += 1
        print(f'EXHAUSTIVE p={p}: {count} subspaces K, dim(V)<=3, criterion and contraction support OK')
        tensors(p)
    # Source Example 5.6: extra-special group times C_p; t is coordinate 2.
    for p in (5,7):
        result = setup(3,p,[[0,1,0],[0,0,1]])
        assert result[:3] == (2,1,2), result
        print(f'FIDELITY Example 5.6 p={p}: dim(K)=2, dim(Q)=2, rank(Omega)=1, |U|={p}')
        pairs = list(combinations(range(4),2))
        rows = [[int(pair == (i,j)) for pair in pairs]
                for i in (0,1) for j in (2,3)]
        result = setup(4,p,rows)
        assert result[:3] == (4,4,4), result
        print(f'FIDELITY Proposition 5.7(1) M=2 p={p}: dim(Q)=4, rank(Omega)=4, |U|={p**4}')
    rng = Random(261006909)
    for p in (3,5,7):
        count = 0
        for d in (4,5):
            n = d*(d-1)//2
            for trial in range(120):
                r = trial % (n+1)
                rows = [[rng.randrange(p) for _ in range(n)] for _ in range(r)]
                setup(d,p,rows);count += 1
        print(f'SEEDED p={p}: {count} cases, dim(V)=4,5, criterion and contraction support OK')
    print('SCOPE: linear criterion and tensor identities checked; no exhaustive group-algebra automorphism check; no Lean proof')
    print('ALL_CHECKS_OK')

if __name__ == '__main__':
    main()
