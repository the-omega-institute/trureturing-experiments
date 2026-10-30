#!/usr/bin/env python3
"""Exact finite evidence for Mizzi 2603.27559v3, Section 7, second assertion.
Python standard library only; no floating point and no universal proof claim.
"""
from collections import Counter
from itertools import combinations, product


def graph(n, edges):
    g = [set() for _ in range(n)]
    for u,v in edges:
        assert 0 <= u < n and 0 <= v < n and u != v
        g[u].add(v)
        g[v].add(u)
    return g


def tf(g,a,b):
    n = len(g)
    assert sorted(a) == sorted(b) == list(range(n))
    return all((v in g[u]) == (b[v] in g[a[u]])
               for u in range(n) for v in range(n))


def inverse(a):
    b = [0]*len(a)
    for u,v in enumerate(a):
        b[v] = u
    return b


def connected(g):
    seen,todo = {0},[0]
    while todo:
        for v in g[todo.pop()]:
            if v not in seen:
                seen.add(v)
                todo.append(v)
    return len(seen) == len(g)


def bipartition(g):
    color = {}
    for root in range(len(g)):
        if root in color:
            continue
        color[root] = 0
        todo = [root]
        while todo:
            u = todo.pop()
            for v in g[u]:
                if v not in color:
                    color[v] = 1-color[u]
                    todo.append(v)
                elif color[v] == color[u]:
                    return None
    return color


def automorphisms(g):
    """Complete enumeration, pruned by exact invariants and adjacency tests."""
    n = len(g)
    colors = [len(x) for x in g]
    while True:
        sigs = [(colors[u],tuple(sorted(colors[v] for v in g[u]))) for u in range(n)]
        labels = {s:i for i,s in enumerate(sorted(set(sigs)))}
        new = [labels[s] for s in sigs]
        stable = len(set(new)) == len(set(colors))
        colors = new
        if stable:
            break
    candidates = [[v for v in range(n) if colors[u] == colors[v]] for u in range(n)]
    mapping,used = {},set()
    def search():
        if len(mapping) == n:
            yield tuple(mapping[u] for u in range(n))
            return
        choices = {u:[v for v in candidates[u] if v not in used and
                      all((x in g[u]) == (y in g[v]) for x,y in mapping.items())]
                   for u in range(n) if u not in mapping}
        u = min(choices,key=lambda x:(len(choices[x]),x))
        for v in choices[u]:
            mapping[u] = v
            used.add(v)
            yield from search()
            used.remove(v)
            del mapping[u]
    return list(search())


def cycles(g):
    """Each unoriented simple cycle once, its least vertex first."""
    result = []
    for root in range(len(g)):
        def dfs(path,seen):
            for v in sorted(g[path[-1]]):
                if v == root:
                    if len(path) >= 3 and path[1] < path[-1]:
                        result.append(tuple(path))
                elif v > root and v not in seen:
                    dfs(path+[v],seen|{v})
        dfs([root],{root})
    return result


def verify_cycle(g,c):
    assert len(c) >= 3 and len(set(c)) == len(c)
    assert all(c[(i+1)%len(c)] in g[c[i]] for i in range(len(c)))


def lifted_components(p,sums):
    k = len(sums)
    h = graph(p*k,((i*p+x,((i+1)%k)*p+(s-x)%p)
                   for i,s in enumerate(sums) for x in range(p)))
    assert all(len(n) == 2 for n in h)
    seen,components = set(),[]
    for root in range(len(h)):
        if root in seen:
            continue
        path,prev,curr = [],None,root
        while curr not in seen:
            path.append(curr)
            seen.add(curr)
            nxt = min(v for v in h[curr] if v != prev)
            prev,curr = curr,nxt
        assert curr == root
        verify_cycle(h,path)
        components.append(path)
    assert sorted(map(len,components)) == [k]+[2*k]*((p-1)//2)
    returns = []
    for x in range(p):
        y = x
        for s in sums:
            y = (s-y)%p
        returns.append(y)
    assert returns == [(returns[0]-x)%p for x in range(p)]
    assert sum(returns[x] == x for x in range(p)) == 1
    return components


def fidelity():
    # TF-Cousins.tex, Figure gammacockroach, including its closed top triangle.
    edges = [(6,7),(7,12),(12,6),(6,1),(6,3),(1,5),(3,4),(5,10),
             (4,9),(9,10),(4,2),(5,2),(4,11),(5,8),(11,8),
             (1,12),(2,11),(3,10)]
    g = graph(12,((u-1,v-1) for u,v in edges))
    a = [3*(u//3)+(u+1)%3 for u in range(12)]
    assert tf(g,a,inverse(a)) and a != inverse(a)
    assert connected(g) and bipartition(g) is None
    assert len({tuple(sorted(n)) for n in g}) == 12
    assert automorphisms(g) == [tuple(range(12))]
    cs = cycles(g)
    for k in [3,5]:
        pairs = [(c,d) for c in cs for d in cs if len(c) == k and
                 len(d) == 2*k and (k != 3 or set(c).isdisjoint(d))]
        assert pairs
        c,d = min(pairs)
        verify_cycle(g,c)
        verify_cycle(g,d)
        print(f'SOURCE12 k={k}: cycles {tuple(x+1 for x in c)} and {tuple(x+1 for x in d)}; disjoint={set(c).isdisjoint(d)}')
    print(f'SOURCE12 automorphisms=1 TF_pair_verified=True cycle_length_counts={dict(sorted(Counter(map(len,cs)).items()))}')
    pin = graph(7,[(i,(i+1)%6) for i in range(6)]+[(6,0),(6,3)])
    assert sorted(set(map(len,cycles(pin)))) == [5,6]
    assert len(automorphisms(pin)) == 4
    print('SOURCE7 cycle_lengths=[5, 6] automorphisms=4')
    old = graph(10,((u-1,v-1) for u,v in [(1,2),(1,10),(2,4),(3,4),
                (3,6),(4,5),(5,7),(5,10),(6,10),(7,9),(8,10)]))
    swap = list(range(10))
    swap[0],swap[5],swap[1],swap[2] = 5,0,2,1
    assert tf(old,swap,swap) and swap != list(range(10))
    assert sorted(map(len,cycles(old))) == [5,5,6]
    print('SRIVASTAVA nonidentity_ordinary_automorphism=(1 6)(2 3) cycle_lengths=[5, 5, 6]')


def small_invariant_graphs():
    total,bip_count,odd_count = 0,0,0
    layouts = [(3,r,f) for r in range(1,4) for f in range(3) if 3*r+f <= 10]
    layouts += [(5,r,f) for r in range(1,3) for f in range(2)]
    for p,r,f in layouts:
        n = p*r+f
        a = [(u//p)*p+(u+1)%p if u < p*r else u for u in range(n)]
        ai = inverse(a)
        unused = set(combinations(range(n),2))
        blocks = []
        while unused:
            todo = [min(unused)]
            block,loop = set(),False
            while todo:
                u,v = todo.pop()
                pair = tuple(sorted((u,v)))
                if u == v:
                    loop = True
                    continue
                if pair in block:
                    continue
                block.add(pair)
                todo.extend([(a[u],ai[v]),(a[v],ai[u])])
            unused -= block
            if not loop:
                blocks.append(block)
        assert len(blocks) <= 12
        for mask in range(1 << len(blocks)):
            g = graph(n,(e for i,b in enumerate(blocks) if mask >> i & 1 for e in b))
            assert tf(g,a,ai)
            assert all(not g[u].intersection(range((u//p)*p,(u//p+1)*p))
                       for u in range(p*r))
            assert all(len(g[v].intersection(range(i*p,(i+1)*p))) in [0,p]
                       for v in range(p*r,n) for i in range(r))
            q = graph(r,((i,j) for i,j in combinations(range(r),2)
                        if any(v in g[u] for u in range(i*p,(i+1)*p)
                               for v in range(j*p,(j+1)*p))))
            color = bipartition(q)
            if color is not None:
                switch = [a[u] if color[u//p] == 0 else ai[u] for u in range(p*r)]
                switch += list(range(p*r,n))
                assert tf(g,switch,switch) and switch != list(range(n))
                bip_count += 1
            else:
                qc = min(c for c in cycles(q) if len(c)%2)
                sums = []
                for i,j in zip(qc,qc[1:]+qc[:1]):
                    x,y = min((x,y) for x in range(p) for y in range(p)
                              if j*p+y in g[i*p+x])
                    sums.append((x+y)%p)
                for component in lifted_components(p,sums):
                    verify_cycle(g,[qc[u//p]*p+u%p for u in component])
                odd_count += 1
            total += 1
    print(f'INVARIANT_GRAPHS layouts={len(layouts)} total={total} bipartite_quotient={bip_count} odd_quotient={odd_count}')


def all_small_matchings():
    cases = 0
    ranges = [(3,3),(3,5),(3,7),(5,3),(5,5),(7,3)]
    for p,k in ranges:
        for sums in product(range(p),repeat=k):
            lifted_components(p,sums)
            cases += 1
    print(f'REFLECTION_LIFTS all_sum_tuples cases={cases} prime_and_length_ranges={ranges}')


def main():
    fidelity()
    small_invariant_graphs()
    all_small_matchings()
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
