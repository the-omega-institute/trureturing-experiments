#!/usr/bin/env python3
"""Exact preregistration checks for Burešová–Pták–Ševic, arXiv:2609.10770v1.

The finite tests are evidence for a proof route, not a proof of the conjecture.
Only integer and Fraction arithmetic enters the tests.
"""
from collections import deque
from fractions import Fraction as F
from itertools import product


def compositions(total, n):
    if n == 1:
        yield (total,)
    else:
        for x in range(total + 1):
            for rest in compositions(total - x, n - 1):
                yield (x,) + rest


def neighbors(edge, rows):
    return {j for i, j in edge if i in rows}


def hall(edge, p, q):
    for mask in range(1 << len(p)):
        rows = {i for i in range(len(p)) if mask >> i & 1}
        if sum(p[i] for i in rows) > sum(q[j] for j in neighbors(edge, rows)):
            return False
    return True


def transport(edge, p, q):
    """Integral augmenting-path flow, independent of the Hall predicate."""
    m, n, total = len(p), len(q), sum(p)
    assert total == sum(q)
    source, sink = m + n, m + n + 1
    cap = [[0] * (m + n + 2) for _ in range(m + n + 2)]
    for i in range(m):
        cap[source][i] = p[i]
    for j in range(n):
        cap[m+j][sink] = q[j]
    for i, j in edge:
        cap[i][m+j] = total
    residual = [row[:] for row in cap]
    value = 0
    while True:
        parent = {source: None}
        queue = deque([source])
        while queue and sink not in parent:
            u = queue.popleft()
            for v in range(len(cap)):
                if residual[u][v] > 0 and v not in parent:
                    parent[v] = u
                    queue.append(v)
        if sink not in parent:
            break
        amount, v = total, sink
        while parent[v] is not None:
            u = parent[v]
            amount = min(amount, residual[u][v])
            v = u
        v = sink
        while parent[v] is not None:
            u = parent[v]
            residual[u][v] -= amount
            residual[v][u] += amount
            v = u
        value += amount
    if value != total:
        return None
    c = [[cap[i][m+j] - residual[i][m+j] for j in range(n)] for i in range(m)]
    assert all(x >= 0 for row in c for x in row)
    assert all(sum(c[i]) == p[i] for i in range(m))
    assert all(sum(c[i][j] for i in range(m)) == q[j] for j in range(n))
    assert all(c[i][j] == 0 for i in range(m) for j in range(n) if (i,j) not in edge)
    return c


def check_transportation():
    total_cases = feasible = 0
    denominator = 4
    for m in range(1, 4):
        for n in range(1, 4):
            pairs = list(product(range(m), range(n)))
            for mask in range(1 << len(pairs)):
                edge = {e for k, e in enumerate(pairs) if mask >> k & 1}
                for p in compositions(denominator, m):
                    for q in compositions(denominator, n):
                        expected = hall(edge, p, q)
                        c = transport(edge, p, q)
                        assert (c is not None) == expected, (m,n,edge,p,q,c)
                        total_cases += 1
                        if c is not None:
                            feasible += 1
                            rational = [[F(x, denominator) for x in row] for row in c]
                            assert sum(map(sum, rational), F(0)) == 1
    print(f'TRANSPORT: {total_cases} graphs/marginals, {feasible} feasible; all m,n=1..3, denominator=4')


def source_blocks():
    """Maximal segments read from simplified_js_2.tex (the source figure)."""
    blocks = []
    for k, letter in enumerate('cdgef'):
        blocks.extend([('a', f'u{k}', f'L{k}'), ('b', f'v{k}', f'R{k}'),
                       (f'L{k}', letter+'1', letter+'2', f'R{k}')])
    for k, (x,y) in enumerate([('c1','d1'), ('g1','e1'), ('c1','f1'),
                              ('c2','e2'), ('c2','g2'), ('d2','f2')]):
        blocks.append((x,y,f'z{k}'))
    return blocks


def check_source():
    blocks = source_blocks()
    atoms = sorted(set().union(*map(set,blocks)))
    state = {a:F(0) for a in atoms}
    state['a'] = state['b'] = F(1)
    for letter in 'cdgef':
        state[letter+'1'] = state[letter+'2'] = F(1,2)
    assert all(0 <= x <= 1 for x in state.values())
    assert all(sum(state[x] for x in block) == 1 for block in blocks)
    assert state['a'] == state['b'] == 1
    # Under a=b=1 the five horizontal equations force x2=1-x1.
    # Each internal third atom must be 1-x-y and nonnegative.
    candidates = 0
    for bits in product((0,1), repeat=5):
        s = dict(state)
        for letter, bit in zip('cdgef', bits):
            s[letter+'1'], s[letter+'2'] = F(bit), F(1-bit)
        for k, block in enumerate(blocks[-6:]):
            s[f'z{k}'] = 1-s[block[0]]-s[block[1]]
        if all(x in (0,1) for x in s.values()) and all(sum(s[x] for x in block)==1 for block in blocks):
            candidates += 1
    assert candidates == 0
    print(f'SOURCE Proposition 2.9: {len(atoms)} atoms, {len(blocks)} blocks; critical half-state verified; no two-valued critical state (32 cases)')
    # The source's two end blocks; m=L0 and n=R0 share a horizontal block.
    A, B = ('a','u0','L0'), ('b','v0','R0')
    orthogonal = {frozenset((x,y)) for block in blocks for x in block for y in block if x != y}
    edge = {(i,j) for i,x in enumerate(A) for j,y in enumerate(B) if frozenset((x,y)) not in orthogonal}
    assert edge == set(product(range(3),repeat=2)) - {(2,2)}
    # The source's forbidden m=n=1 pair fails the cut at {m}.
    assert not hall(edge, (0,0,1), (0,0,1))
    assert transport(edge, (0,0,1), (0,0,1)) is None
    p, q = (2,1,1), (2,1,1)
    c = transport(edge,p,q)
    assert hall(edge,p,q) and c is not None
    print('SOURCE end-block transport at marginals (1/2,1/4,1/4):', [[str(F(x,4)) for x in row] for row in c])


def partitions(n):
    def rec(k, labels, largest):
        if k == n:
            yield tuple(frozenset(i for i,x in enumerate(labels) if x==j) for j in range(largest+1))
            return
        for j in range(largest+2):
            yield from rec(k+1,labels+[j],max(largest,j))
    if n:
        yield from rec(1,[0],0)


def check_cut_reduction():
    """All two Boolean subalgebras of powersets on <=4 points, all atom subsets.

    A fiber is an atom block; ambient order is actual subset inclusion.
    Nonempty atom intersections supply Dirac global states on allowed pairs.
    """
    cut_cases = marginal_cases = 0
    for n in range(1,5):
        universe = frozenset(range(n))
        ps = list(partitions(n))
        for A in ps:
            for B in ps:
                edge = {(i,j) for i,a in enumerate(A) for j,b in enumerate(B) if a & b}
                for mask in range(1 << len(A)):
                    rows = {i for i in range(len(A)) if mask >> i & 1}
                    x = frozenset().union(*(A[i] for i in rows))
                    ns = neighbors(edge,rows)
                    y = frozenset().union(*(B[j] for j in ns))
                    assert x <= y
                    assert all(B[j] <= universe-x for j in range(len(B)) if j not in ns)
                    cut_cases += 1
                for weights in compositions(4,n):
                    p = tuple(sum(weights[x] for x in a) for a in A)
                    q = tuple(sum(weights[x] for x in b) for b in B)
                    c = transport(edge,p,q)
                    assert hall(edge,p,q) and c is not None
                    # Lift the coupling to Dirac states, checking every event in both algebras.
                    lifted = [F(0) for _ in range(n)]
                    for i,j in edge:
                        lifted[min(A[i] & B[j])] += F(c[i][j],4)
                    assert sum(lifted) == 1
                    for block,mass in list(zip(A,p))+list(zip(B,q)):
                        assert sum(lifted[x] for x in block) == F(mass,4)
                    marginal_cases += 1
    print(f'CUT/LIFT: {cut_cases} subset cuts and {marginal_cases} rational common-source marginals; all partition pairs on 1..4 points')


def main():
    check_source()
    check_transportation()
    check_cut_reduction()
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
